function F = emg_features(x, fs, varargin)
%EMG_FEATURES  Windowed feature extraction for sEMG classification.
%
%   F = EMG_FEATURES(x, fs)
%   F = EMG_FEATURES(x, fs, 'win_ms', 200, 'inc_ms', 50, 'labels', L)
%
%   x   [nsamples x nchannels]  FILTERED EMG (run emg_filter first)
%
%   OUTPUT
%     F.X        [nwin x (nfeat*nchannels)]  feature matrix
%     F.names    feature names, in column order
%     F.t        window centre times, s
%     F.y        majority label per window, if 'labels' was supplied
%     F.win, F.inc   the window and increment actually used, in samples
%
%   ---------------------------------------------------------------------------
%   FEATURE SET, AND WHY EACH ONE IS HERE
%
%   LITERATURE-VALIDATED 2026-09-25 - see DECISION_LOG D-25 for the full pass.
%
%   TIME DOMAIN - the Hudgins set, the standard baseline in myoelectric control.
%     Hudgins B, Parker P, Scott RN (1993). "A new strategy for multifunction
%     myoelectric control." IEEE Trans. Biomed. Eng. 40(1):82-94.
%     doi:10.1109/10.204774. RMS/VAR are the standard amplitude add-ons to that
%     original four-feature set, not part of Hudgins' own four.
%     RMS   root mean square          amplitude; correlates with force
%     MAV   mean absolute value       amplitude, cheaper than RMS (Hudgins)
%     WL    waveform length           cumulative path length: amplitude AND
%                                     frequency content in one number (Hudgins)
%     ZC    zero crossings            frequency proxy, threshold-gated (Hudgins)
%     SSC   slope sign changes        frequency proxy, threshold-gated (Hudgins)
%     VAR   variance                  amplitude
%
%   FREQUENCY DOMAIN - and these matter specifically for THIS project
%     MNF   mean frequency
%     MDF   median frequency
%
%   Post-stroke sEMG shows a documented DOWNWARD SPECTRAL SHIFT (loss of
%   high-threshold fast-twitch motor units). MNF and MDF capture that directly.
%   Citation with numbers: Yao B, Zhang X, Li S, Li X, Chen X, Klein CS, Zhou P
%   (2015). "Analysis of linear electrode array EMG for assessment of
%   hemiparetic biceps brachii muscles." Front. Hum. Neurosci. 9:569.
%   Reported paretic vs contralateral: median frequency 64.56+/-1.32 Hz vs
%   79.46+/-2.45 Hz (p<0.001); mean power frequency 71.36+/-1.58 Hz vs
%   86.60+/-1.67 Hz (p<0.001) - both significantly lower on the paretic side.
%   (That study is on biceps brachii, not the forearm muscles used here; the
%   mechanism - selective loss of high-threshold fast-twitch motor units is
%   general, but the exact numbers shall not be quoted as forearm-specific in
%   the thesis without a forearm-muscle source, if a more specific one is found.)
%   An amplitude-only feature set cannot distinguish "weak because few units are
%   firing" from "weak because the electrode is badly placed"; adding the
%   spectral features gives the classifier something that amplitude alone
%   does not carry.  Keep them even if they cost accuracy on healthy subjects.
%
%   WINDOWING.  200 ms window, 50 ms increment by default.  Smith LH, Hargrove
%   LJ, Lock BA, Kuiken TA (2011), "Determining the optimal window length for
%   pattern recognition-based myoelectric control: balancing the competing
%   effects of classification error and controller delay," IEEE Trans. Neural
%   Syst. Rehabil. Eng., tested 50-550 ms windows and found 150-250 ms to be
%   the sweet spot balancing classification error against controller delay.
%   This project's 200 ms sits inside that band, not just under a generic
%   "300 ms feels responsive" ceiling (that ceiling itself traces to Farrell &
%   Weir (2007), "The Optimal Controller Delay for Myoelectric Prostheses,"
%   IEEE Trans. Neural Syst. Rehabil. Eng., and Englehart & Hudgins (2003),
%   "A Robust, Real-Time Control Scheme for Multifunction Myoelectric
%   Control," IEEE Trans. Biomed. Eng. 50(7):848-854).  The 50 ms increment
%   gives a 20 Hz decision rate, comfortably above the 50 Hz intent contract
%   after zero-order hold.
%
%   ZC and SSC use a DEADZONE threshold. Without it they count noise crossings
%   and the feature becomes a noise meter.  Default is 3x the per-channel
%   resting RMS, adaptive per session rather than a fixed universal value; the
%   literature on this specifically (e.g. the multi-dataset threshold-selection
%   studies around Phinyomark-era EMG feature work) generally warns that a
%   single FIXED threshold across subjects/sessions degrades performance,
%   which is the argument for keeping this adaptive rather than hard-coding a
%   number -to supply 'thresh' explicitly for reproducibility across sessions.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev B 2026-09-25

o.win_ms = 200;
o.inc_ms = 50;
o.labels = [];
o.thresh = [];          % ZC/SSC deadzone per channel; [] -> 3x resting RMS
o.rest_win = [];
o.channels = {};
o = parseopts(o, varargin);

[N, nch] = size(x);
if isempty(o.channels)
    o.channels = arrayfun(@(k)sprintf('ch%d',k), 1:nch, 'uni', 0);
end

win = max(8, round(o.win_ms*fs/1000));
inc = max(1, round(o.inc_ms*fs/1000));
starts = 1:inc:(N-win+1);
nw = numel(starts);
if nw < 1
    error('emg_features:tooShort','Signal is shorter than one window.');
end

% ---- ZC / SSC deadzone ------------------------------------------------------
if isempty(o.thresh)
    if isempty(o.rest_win)
        rr = 1:min(round(0.5*fs), N);
    else
        rr = o.rest_win(1):o.rest_win(2);
    end
    o.thresh = 3 * sqrt(mean(x(rr,:).^2, 1));
end
o.thresh = o.thresh(:)';
if numel(o.thresh) == 1, o.thresh = repmat(o.thresh, 1, nch); end

featnames = {'RMS','MAV','WL','ZC','SSC','VAR','MNF','MDF'};
nf = numel(featnames);
X  = zeros(nw, nf*nch);

nfft = 2^nextpow2(win);
f    = (0:nfft/2)' * fs/nfft;
w    = 0.5 - 0.5*cos(2*pi*(0:win-1)'/(win-1));      % Hann

for k = 1:nw
    seg = x(starts(k):starts(k)+win-1, :);
    for c = 1:nch
        v  = seg(:,c);
        dv = diff(v);
        th = o.thresh(c);

        rmsv = sqrt(mean(v.^2));
        mav  = mean(abs(v));
        wl   = sum(abs(dv));
        zc   = sum( (v(1:end-1).*v(2:end) < 0) & (abs(dv) > th) );
        ssc  = sum( (dv(1:end-1).*dv(2:end) < 0) & ...
                    (max(abs(dv(1:end-1)),abs(dv(2:end))) > th) );
        varv = var(v);

        P = abs(fft((v-mean(v)).*w, nfft)).^2;
        P = P(1:nfft/2+1);
        P(2:end-1) = 2*P(2:end-1);
        tot = sum(P);
        if tot > 0
            mnf = sum(f.*P)/tot;
            cum = cumsum(P);
            mdf = f(find(cum >= 0.5*tot, 1, 'first'));
        else
            mnf = 0; mdf = 0;
        end

        X(k, (c-1)*nf + (1:nf)) = [rmsv mav wl zc ssc varv mnf mdf];
    end
end

names = cell(1, nf*nch);
for c = 1:nch
    for j = 1:nf
        names{(c-1)*nf + j} = sprintf('%s_%s', o.channels{c}, featnames{j});
    end
end

F.X = X;  F.names = names;  F.featnames = featnames;
F.t = (starts + win/2 - 1)'/fs;
F.win = win;  F.inc = inc;  F.thresh = o.thresh;
F.channels = o.channels;

% ---- majority label per window ---------------------------------------------
if ~isempty(o.labels)
    L = o.labels(:);
    y = zeros(nw,1);
    for k = 1:nw
        y(k) = mode(L(starts(k):starts(k)+win-1));
    end
    F.y = y;
end

if nargout == 0 || true
    fprintf('\n  EMG FEATURES\n');
    fprintf('  ------------------------------------------------------\n');
    fprintf('  window %d ms (%d samples), increment %d ms -> %.0f Hz\n', ...
            o.win_ms, win, o.inc_ms, 1000/o.inc_ms);
    fprintf('  %d windows x %d features (%d per channel x %d channels)\n', ...
            nw, size(X,2), nf, nch);
    fprintf('  features: %s\n', strjoin(featnames, ', '));
    fprintf('  ZC/SSC deadzone: %s\n', ...
            strjoin(arrayfun(@(v)sprintf('%.4g',v), o.thresh, 'uni',0), ' '));
    if isfield(F,'y')
        u = unique(F.y);
        fprintf('  labelled: %d classes, counts', numel(u));
        for k = 1:numel(u), fprintf(' %d:%d', u(k), sum(F.y==u(k))); end
        fprintf('\n');
    end
    fprintf('  ------------------------------------------------------\n');
    fprintf('  MNF/MDF are in the set deliberately - post-stroke EMG shifts\n');
    fprintf('  DOWN in frequency, and amplitude features cannot see that.\n\n');
end
end

% ===============================================================================
function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_features:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
