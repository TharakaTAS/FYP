function Q = emg_quality(x, fs, varargin)
%EMG_QUALITY  Signal-quality gate for sEMG.  RUN THIS BEFORE ANY PROCESSING.
%
%   Q = EMG_QUALITY(x, fs)
%   Q = EMG_QUALITY(x, fs, 'rest_win', [i1 i2], 'channels', {'FCR','ECR'}, ...
%                   'mains', 50, 'plot', true)
%
%   x   [nsamples x nchannels] RAW EMG, before any filtering.
%
%   ---------------------------------------------------------------------------
%   WHY THIS EXISTS, AND WHY IT COMES FIRST
%
%   This is the stage students skip and reviewers attack.  A beautifully filtered
%   signal from a loose electrode is still garbage, and no amount of downstream
%   cleverness recovers it.  Worse: bad channels fail SILENTLY - you get numbers,
%   they just do not mean anything.
%
%   The rule for this project: run this at CAPTURE TIME, while the subject is
%   still in the chair.  A failed channel is then a two-minute re-prep.  A week
%   later it is an unusable session.
%
%   CHECKS, and what a failure actually means physically
%     clip_pct      samples at the ADC rail        -> gain too high, or artefact
%     flat_pct      near-zero variance runs        -> disconnected lead
%     snr_db        active RMS vs resting RMS      -> poor skin prep, or truly
%                                                     no voluntary activity
%     mains_ratio   power at 50 Hz +/- 1 vs total  -> poor grounding, bad contact
%     drift_pct     sub-2 Hz power vs total        -> motion artefact, cable sway
%     dc_offset     mean level                     -> half-cell potential drift
%     rest_rms      the resting noise floor        -> sets the intent deadband
%
%   OUTPUT Q has one row per channel, a .pass flag per channel, and .allpass.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev A 2026-09-20

o.rest_win = [];          % [i1 i2] indices of a quiet rest period
o.channels = {};
o.mains    = 50;
o.vrange   = 5;           % volts, ADC full scale (+/-), for clipping detection
o.plot     = false;
% thresholds
o.th_clip  = 0.1;         % % of samples
o.th_flat  = 1.0;         % %
o.th_snr   = 10;          % dB
o.th_mains = 0.05;        % fraction of total power
o.th_drift = 0.60;        % low-band power as a fraction of low + EMG band
o.mode     = 'live';      % 'live' = capture-time gate (strict)
                          % 'archival' = someone else's recorded file: drift is
                          %   reported but not failed, because you cannot go
                          %   back and re-tape their cable, and the band-pass
                          %   removes it anyway.
o = parseopts(o, varargin);

[N, nch] = size(x);
if isempty(o.channels)
    o.channels = arrayfun(@(k)sprintf('ch%d',k), 1:nch, 'uni', 0);
end
if isempty(o.rest_win)
    % No rest window given - use the quietest 1 s as a proxy, and say so.
    win = min(round(fs), N);
    e   = movmean(mean(x.^2,2), win);
    [~, c] = min(e(win:end-win));
    o.rest_win = [c, c+win-1];
    restNote = 'auto (quietest 1 s)';
else
    restNote = 'supplied';
end
r = o.rest_win(1):o.rest_win(2);

Q = struct();
Q.channels = o.channels(:);
Q.rest_source = restNote;
Q.pass = false(nch,1);

for c = 1:nch
    v = x(:,c);

    % --- clipping
    % Clipping means samples PINNED at a rail, i.e. many samples sharing the
    % exact same extreme value.  A large DC OFFSET is not clipping - the old
    % test (|v| >= 0.995*vrange) reported 100% clipping on a channel that was
    % merely sitting at -5 V, which is a half-cell offset, not saturation.
    vmax = max(v);  vmin = min(v);
    tol  = 1e-9 + 1e-6*(vmax - vmin);
    clip_pct = 100 * mean( abs(v-vmax) < tol | abs(v-vmin) < tol );

    % --- flatline (run of near-zero variance)
    mv = movstd(v, round(0.05*fs));
    flat_pct = 100 * mean(mv < 1e-6);

    % --- SNR: active vs rest
    % The old version compared the supplied rest window against EVERYTHING
    % else - but in a 640 s record with 50 short contractions, "everything
    % else" is mostly rest too, so the ratio came out near 1 and every channel
    % failed.  Instead: band-pass, take an envelope, and compare its top
    % decile (contraction) against its bottom decile (true rest).  That is
    % self-calibrating and does not need to be told where the rest is.
    [bq, aq] = butter(2, [20 min(450, 0.45*fs)]/(fs/2), 'bandpass');
    vb  = filtfilt(bq, aq, v);
    ev  = sqrt(movmean(vb.^2, max(3, round(0.2*fs))));
    evs = sort(ev);
    n10 = max(1, round(0.10*numel(evs)));
    rest_rms = mean(evs(1:n10));
    act_rms  = mean(evs(end-n10+1:end));
    snr_db   = 20*log10(max(act_rms,eps) / max(rest_rms,eps));

    % --- spectrum-based checks
    nfft = 2^nextpow2(min(N, 4*fs));
    [pxx, f] = pwelch_local(v, nfft, fs);
    tot = trapz(f, pxx);
    if ~isempty(o.mains)
        bm = f >= o.mains-1 & f <= o.mains+1;
        mains_ratio = trapz(f(bm), pxx(bm)) / max(tot,eps);
    else
        mains_ratio = 0;
    end
    % Drift: power BELOW the EMG band relative to power INSIDE it.  Measuring
    % it against TOTAL power is meaningless, because the sub-20 Hz baseline
    % wander usually IS most of the total - which is exactly why we band-pass.
    % What matters is whether it swamps the signal we intend to keep.
    bd  = f <= 20;
    bemg = f > 20 & f <= 450;
    p_lo  = trapz(f(bd),  pxx(bd));
    p_emg = trapz(f(bemg), pxx(bemg));
    drift_pct = p_lo / max(p_lo + p_emg, eps);

    dc_offset = mean(v);

    % --- verdict
    fails = {};
    if clip_pct    > o.th_clip,  fails{end+1} = 'CLIPPING';      end %#ok<AGROW>
    if flat_pct    > o.th_flat,  fails{end+1} = 'FLATLINE';      end %#ok<AGROW>
    if snr_db      < o.th_snr,   fails{end+1} = 'LOW SNR';       end %#ok<AGROW>
    if mains_ratio > o.th_mains, fails{end+1} = 'MAINS';         end %#ok<AGROW>
    if drift_pct   > o.th_drift && strcmpi(o.mode,'live')
        fails{end+1} = 'DRIFT';                                   %#ok<AGROW>
    end

    Q.clip_pct(c,1)    = clip_pct;
    Q.flat_pct(c,1)    = flat_pct;
    Q.snr_db(c,1)      = snr_db;
    Q.mains_ratio(c,1) = mains_ratio;
    Q.drift_pct(c,1)   = drift_pct;
    Q.dc_offset(c,1)   = dc_offset;
    Q.rest_rms(c,1)    = rest_rms;
    Q.fails{c,1}       = fails;
    Q.pass(c,1)        = isempty(fails);
end
Q.allpass = all(Q.pass);

% ---- report ----------------------------------------------------------------
fprintf('\n  EMG SIGNAL QUALITY  (mode: %s, rest: %s)\n', o.mode, restNote);
fprintf('  ---------------------------------------------------------------------\n');
fprintf('  %-6s %7s %7s %8s %7s %7s   %s\n', ...
        'ch','clip%','flat%','SNR dB','mains','drift','verdict');
for c = 1:nch
    if Q.pass(c)
        v = 'PASS';
    else
        % NOTE: do NOT build this with [ ... ] concatenation.  A quote directly
        % after a comma inside square brackets is parsed as a transpose operator
        % and the file will not even compile.  sprintf avoids the ambiguity.
        v = sprintf('FAIL: %s', strjoin(Q.fails{c}, ', '));
    end
    fprintf('  %-6s %7.2f %7.2f %8.1f %7.3f %7.3f   %s\n', ...
            o.channels{c}, Q.clip_pct(c), Q.flat_pct(c), Q.snr_db(c), ...
            Q.mains_ratio(c), Q.drift_pct(c), v);
end
fprintf('  ---------------------------------------------------------------------\n');
if Q.allpass
    fprintf('  ALL CHANNELS PASS.  Proceed.\n');
else
    fprintf('  >> DO NOT PROCEED.  Re-prep the failing channels NOW, while the\n');
    fprintf('     subject is still in the chair.  Remedies:\n');
    fprintf('       LOW SNR / MAINS -> abrade skin (NuPrep), re-wipe, re-gel,\n');
    fprintf('                          check the reference electrode on bone\n');
    fprintf('       DRIFT           -> secure the cable, strain-relief the lead\n');
    fprintf('       CLIPPING        -> reduce amplifier gain\n');
    fprintf('       FLATLINE        -> the lead is off\n');
end
fprintf('  Resting RMS per channel sets the intent deadband (see emg_baseline).\n\n');

if o.plot
    figure('Name','EMG signal quality','Color','w');
    for c = 1:nch
        subplot(nch,1,c);
        plot((0:N-1)/fs, x(:,c), 'LineWidth',0.4); grid on; box on;
        hold on; xline(r(1)/fs,':'); xline(r(end)/fs,':');
        ylabel(o.channels{c});
        if c == 1, title('Raw EMG - dotted lines mark the rest window'); end
        if c == nch, xlabel('time [s]'); end
    end
end
end

% =============================================================================
function [pxx, f] = pwelch_local(v, nfft, fs)
% Welch PSD without the Signal Processing Toolbox.
v = v(:) - mean(v);
w = hann_local(nfft);
step = floor(nfft/2);
nseg = max(1, floor((numel(v)-nfft)/step) + 1);
pxx = zeros(nfft/2+1, 1);
for k = 1:nseg
    i0 = (k-1)*step + 1;
    s  = v(i0:i0+nfft-1) .* w;
    P  = abs(fft(s)).^2;
    pxx = pxx + P(1:nfft/2+1);
end
U   = sum(w.^2);
pxx = pxx / (nseg * fs * U);
pxx(2:end-1) = 2*pxx(2:end-1);
f = (0:nfft/2)' * fs/nfft;
end

function w = hann_local(n)
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
end

function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_quality:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
