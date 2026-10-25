function D = emg_synth(varargin)
%EMG_SYNTH  Synthesis of multi-channel sEMG with a post-stroke ("paretic") mode.
%
%   D = EMG_SYNTH()
%   D = EMG_SYNTH('fs',2000,'condition','paretic','nrep',6)
%
%   REPLACES the lost pathological-EMG emulator, and lets the whole pipeline be
%   developed and validated TODAY with no hardware and no dataset.  Every signal
%   it produces has a KNOWN label, a known onset time and a known amplitude, so
%   any downstream stage can be checked against ground truth rather than against
%   an opinion.
%
%   OPTIONS
%     'fs'        sample rate, Hz                      default 2000
%     'channels'  cellstr of muscle names              default {'FCR','FCU','ECR','ECU'}
%     'motions'   cellstr of motion classes            default {'rest','flex','ext','rdev','udev'}
%     'nrep'      repetitions of each motion           default 5
%     'tmove'     seconds of contraction per rep       default 3
%     'trest'     seconds of rest between reps         default 3
%     'condition' 'healthy' | 'paretic' | 'flaccid'    default 'healthy'
%     'severity'  0..1, how affected (paretic only)    default 0.6
%     'mains'     mains frequency, Hz ([] = none)      default 50
%     'seed'      rng seed                             default 4206
%
%   OUTPUT STRUCT D
%     .x          [nsamples x nchannels] raw EMG, volts
%     .t          time vector, s
%     .fs         sample rate
%     .labels     [nsamples x 1] motion class index per sample  <- GROUND TRUTH
%     .classes    the class names
%     .channels   channel names
%     .onsets     [nrep*nmotion x 2] start/stop sample of each contraction
%     .truth      struct of the generative parameters actually used
%
%   ---------------------------------------------------------------------------
%   WHAT "PARETIC" MEANS HERE, AND WHAT IT DOES NOT
%
%   Post-stroke sEMG differs from healthy sEMG in three documented ways, and this
%   function models all three:
%     1. AMPLITUDE is reduced       - fewer motor units recruited
%     2. The SPECTRUM shifts DOWN   - loss of high-threshold, fast-twitch units
%     3. CO-ACTIVATION increases    - agonist and antagonist fire together
%
%   LITERATURE CHECK 2026-09-25 (DECISION_LOG D-25).  (1) and (2) are supported
%   by Yao et al. (2015), Front. Hum. Neurosci. 9:569 (paretic vs contralateral
%   median/mean-power frequency both significantly lower, biceps brachii; see
%   emg_features.m header for the numbers) and the general fewer-motor-unit
%   mechanism it describes.  (3) is GENUINELY MIXED in the literature and
%   should not be stated as settled in the thesis:
%     - Fellows, Kaus, Ross & Thilmann (1994), Electroencephalogr. Clin.
%       Neurophysiol. 93(2):106-112, found NO exaggerated antagonist
%       co-activation in spastic hemiparesis under isometric elbow torque -
%       "hemiparesis...cannot, under isometric conditions, be attributed to
%       hyperactivity of antagonist muscles."
%     - Li et al. (2022), Front. Neurosci. (co-contraction/corticospinal-tract
%       paper), report a significantly higher co-contraction index on the
%       affected side (0.69+/-0.13) vs unaffected (0.46+/-0.18) and healthy
%       (0.44+/-0.11), correlated with corticospinal-tract impairment.
%   Reconciliation offered here, not claimed as proven: increased co-activation
%   ("spastic co-contraction") is well documented as a real phenomenon in the
%   broader spasticity literature, but its presence/magnitude appears to
%   depend on joint, task (isometric vs dynamic) and measurement window, which
%   is consistent with one older, narrower isometric-elbow study finding none.
%   The synthetic model's monotonic coact-rises-with-severity assumption is a
%   simplification, defensible for pipeline testing, but should NOT be quoted
%   in the thesis as an uncontested finding: cite both sides if this
%   assumption is discussed there.
%
%   It is a SIMULATION for developing and testing the pipeline.  It is NOT
%   patient data and no result from it may be reported as a clinical finding.
%   Said explicitly wherever it appears in the thesis.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev B 2026-09-25

% ---- options ---------------------------------------------------------------
o.fs        = 2000;
o.channels  = {'FCR','FCU','ECR','ECU'};
o.motions   = {'rest','flex','ext','rdev','udev'};
o.nrep      = 5;
o.tmove     = 3;
o.trest     = 3;
o.condition = 'healthy';
o.severity  = 0.6;
o.mains     = 50;
o.seed      = 4206;
o = parseopts(o, varargin);
rng(o.seed);

nch = numel(o.channels);
nmo = numel(o.motions);

% ---- muscle activation matrix ----------------------------------------------
% Rows = motion, columns = channel.  Which muscles fire for which movement.
% Every wrist motion is a muscle PAIR - this is why two EMG channels cannot
% separate all four motions (D-13).
A = zeros(nmo, nch);
idx = @(nm) find(strcmpi(o.channels, nm));
set_ = @(M,ch,v) deal(ch, v);
for k = 1:nmo
    switch lower(o.motions{k})
        case 'rest'                                        % nothing
        case 'flex', A(k, [idx('FCR') idx('FCU')]) = [1.0 0.9];
        case 'ext',  A(k, [idx('ECR') idx('ECU')]) = [1.0 0.9];
        case 'rdev', A(k, [idx('FCR') idx('ECR')]) = [0.8 0.9];
        case 'udev', A(k, [idx('FCU') idx('ECU')]) = [0.9 0.8];
        otherwise,   A(k,:) = 0.5;                         % unknown -> generic
    end
end
A = A(:, 1:nch);
A(isnan(A)) = 0;

% ---- condition-dependent modifiers -----------------------------------------
s = max(0, min(1, o.severity));
switch lower(o.condition)
    case 'healthy'
        amp_gain   = 1.00;   % amplitude scaling
        f_shift    = 1.00;   % spectral centre multiplier
        coact      = 0.05;   % antagonist co-activation fraction
        var_rep    = 0.08;   % rep-to-rep amplitude variability
    case 'paretic'
        amp_gain   = 1 - 0.65*s;      % 35% of healthy at s = 1
        f_shift    = 1 - 0.35*s;      % spectrum shifts down
        coact      = 0.05 + 0.45*s;   % co-activation rises sharply
        var_rep    = 0.08 + 0.25*s;   % much more variable
    case 'flaccid'
        amp_gain   = 0.06;            % barely above noise
        f_shift    = 0.55;
        coact      = 0.10;
        var_rep    = 0.50;
    otherwise
        error('emg_synth:badCondition','condition must be healthy|paretic|flaccid');
end

% ---- build the label timeline ----------------------------------------------
seg = [];  labels_seg = [];
order = repmat(2:nmo, 1, o.nrep);          % skip class 1 (rest) as a "movement"
order = order(randperm(numel(order)));
for k = 1:numel(order)
    seg(end+1)        = o.trest;  labels_seg(end+1) = 1;        %#ok<AGROW> rest
    seg(end+1)        = o.tmove;  labels_seg(end+1) = order(k); %#ok<AGROW>
end
seg(end+1) = o.trest;  labels_seg(end+1) = 1;

nseg   = round(seg * o.fs);
N      = sum(nseg);
t      = (0:N-1)'/o.fs;
labels = zeros(N,1);
onsets = zeros(numel(order), 2);
p = 1;  oi = 0;
for k = 1:numel(nseg)
    r = p:(p+nseg(k)-1);
    labels(r) = labels_seg(k);
    if labels_seg(k) > 1
        oi = oi + 1;  onsets(oi,:) = [r(1) r(end)];
    end
    p = r(end) + 1;
end

% ---- per-channel drive envelope --------------------------------------------
drive = zeros(N, nch);
rampn = round(0.25*o.fs);                       % 250 ms recruitment ramp
w = [linspace(0,1,rampn) ones(1,1)];            %#ok<NASGU>
for k = 1:numel(order)
    r  = onsets(k,1):onsets(k,2);
    cls = order(k);
    g   = amp_gain * (1 + var_rep*randn);       % rep-to-rep variability
    g   = max(g, 0.02);
    env = trapzwin(numel(r), rampn);
    for c = 1:nch
        a = A(cls,c);
        if a == 0, a = coact; end               % antagonist co-activation
        drive(r,c) = drive(r,c) + g*a*env(:);
    end
end

% ---- synthesise the raw signal ---------------------------------------------
% sEMG is well approximated as amplitude-modulated, band-shaped noise.  The
% band centre moves with f_shift to model the post-stroke spectral shift.
x = zeros(N, nch);
base_lo = 20;  base_hi = 450;
hi = min(base_hi*f_shift, 0.45*o.fs);
lo = max(base_lo*f_shift, 5);
[bb, aa] = butter(4, [lo hi]/(o.fs/2), 'bandpass');
for c = 1:nch
    carrier = filtfilt(bb, aa, randn(N,1));
    carrier = carrier / std(carrier);
    baseline = 0.008 * randn(N,1);                        % instrumentation noise
    x(:,c) = (0.55*drive(:,c)) .* carrier + baseline;
end

% ---- artefacts -------------------------------------------------------------
if ~isempty(o.mains)
    for c = 1:nch
        ph = 2*pi*rand;
        x(:,c) = x(:,c) + 0.006*sin(2*pi*o.mains*t + ph) ...
                        + 0.002*sin(2*pi*3*o.mains*t + ph);   % 3rd harmonic
    end
end
% slow motion artefact / baseline wander
[bw, aw] = butter(2, 2/(o.fs/2), 'low');
for c = 1:nch
    x(:,c) = x(:,c) + 0.01*filtfilt(bw, aw, randn(N,1));
end

% ---- output ----------------------------------------------------------------
D.x = x;  D.t = t;  D.fs = o.fs;
D.labels = labels;  D.classes = o.motions;  D.channels = o.channels;
D.onsets = onsets;  D.order = order(:);
D.truth = struct('A',A,'amp_gain',amp_gain,'f_shift',f_shift, ...
                 'coact',coact,'var_rep',var_rep,'opt',o);

if nargout == 0
    fprintf('\n  emg_synth  -  %s (severity %.2f)\n', o.condition, s);
    fprintf('  ------------------------------------------------------\n');
    fprintf('  %d channels: %s\n', nch, strjoin(o.channels,', '));
    fprintf('  %d classes:  %s\n', nmo, strjoin(o.motions,', '));
    fprintf('  %d contractions, %.1f s total at %g Hz\n', ...
            numel(order), N/o.fs, o.fs);
    fprintf('  amplitude x%.2f, spectrum x%.2f, co-activation %.0f%%\n', ...
            amp_gain, f_shift, 100*coact);
    fprintf('  ------------------------------------------------------\n');
    fprintf('  SIMULATION ONLY.  Not patient data.  No result from this\n');
    fprintf('  may be reported as a clinical finding.\n\n');
    clear D
end
end

% =============================================================================
function w = trapzwin(n, rampn)
rampn = min(rampn, floor(n/2));
w = ones(n,1);
if rampn > 1
    w(1:rampn)        = linspace(0,1,rampn);
    w(end-rampn+1:end) = linspace(1,0,rampn);
end
% physiological tremor / firing-rate fluctuation
w = w .* (1 + 0.12*sin(2*pi*(0:n-1)'/n*3.7) + 0.06*randn(n,1));
w = max(w, 0);
end

function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_synth:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
