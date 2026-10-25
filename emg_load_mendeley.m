function D = emg_load_mendeley(csvfile, varargin)
%EMG_LOAD_MENDELEY  Load + segment + IDENTIFY CHANNELS for the Mendeley
%                   ckwc76xr2z 4-channel sEMG hand-gesture dataset.
%
%   D = EMG_LOAD_MENDELEY('10_raw.csv')
%   D = EMG_LOAD_MENDELEY(f, 'fs',2000, 'plot',true)
%
%   DATASET: "Dataset for multi-channel surface electromyography (sEMG) signals
%   of hand gestures", Mendeley Data doi:10.17632/ckwc76xr2z.2, described in
%   Data in Brief (2022), ScienceDirect S2352340922001330.
%
%   WHY THIS DATASET IS A GOOD MATCH FOR MC 4206
%     - Recorded on a BIOPAC MP36 - the same instrument we own
%     - Ag/AgCl bipolar surface electrodes - our electrode type
%     - The four muscles are ECR, FCR, ECU, FCU - exactly our montage (D-13)
%     - It contains all four of our wrist motions: extension, flexion, ulnar
%       deviation, radial deviation (plus six we do not use)
%     - 2 kHz, 640 s per participant, 40 participants
%
%   PUBLISHED PROTOCOL (from the Data in Brief paper, verified against the files)
%     6 s gesture hold, 4 s rest before and after, 30 s rest between cycles,
%     5 cycles x 10 gestures, 640 s total, fs = 2 kHz.
%     Gesture order: rest, extension, flexion, ulnar dev, radial dev, grip,
%                    abduction, adduction, supination, pronation.
%
%   ---------------------------------------------------------------------------
%   THE PROBLEM THIS FUNCTION SOLVES
%
%   The raw CSV is four unlabelled numeric columns.  No header, no timestamp,
%   no marker channel, no label column.  Neither the Mendeley page nor the paper
%   states WHICH COLUMN IS WHICH MUSCLE.
%
%   The advice circulating online - "just take the FCR and ECR columns" - cannot
%   be followed as written, because nothing tells you which columns those are.
%   Guessing has a 1-in-12 chance of getting the pair right, and a wrong guess
%   does not announce itself: you still get a classifier, still get an accuracy,
%   and the number is meaningless.
%
%   So this function identifies the channels from the data, using physiology:
%
%       FCR  active in FLEXION   and RADIAL deviation
%       FCU  active in FLEXION   and ULNAR  deviation
%       ECR  active in EXTENSION and RADIAL deviation
%       ECU  active in EXTENSION and ULNAR  deviation
%
%   Those four gestures give each muscle a unique two-bit signature.  We score
%   all 24 possible column-to-muscle assignments against the measured activity
%   and report the best, with a margin over the runner-up.  A small margin means
%   DO NOT TRUST THE MAPPING - say so rather than proceeding.
%
%   This is a stronger thing to put in a thesis than a forum post, and it costs
%   nothing: the data identifies itself.
%
%   OUTPUT D
%     .x        [N x 4] raw, volts          .fs, .t
%     .env      [N x 4] RMS envelope        .Xf  filtered
%     .seg      [nseg x 3]  start, stop, gesture index
%     .gestures gesture names
%     .map      1x4 struct: which column is which muscle, with confidence
%     .Q        signal-quality report (emg_quality)
%     .profile  [10 x 4] mean normalised activity, gesture x channel
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev A 2026-09-20

% Make this script find its own neighbours (emg_quality, emg_filter, ...)
% no matter which folder MATLAB happens to be sitting in.
addpath(fileparts(mfilename('fullpath')));

o.fs        = 2000;
o.t_hold    = 6;      % s, gesture
o.t_rest    = 4;      % s, rest each side
o.t_cycle   = 30;     % s, between cycles
o.ncycle    = 5;
o.plot      = true;
o.gestures  = {'rest','extension','flexion','ulnar_dev','radial_dev', ...
               'grip','abduction','adduction','supination','pronation'};
o = parseopts(o, varargin);
fs = o.fs;

% ---- load ------------------------------------------------------------------
x = readmatrix(csvfile);
if size(x,2) ~= 4
    error('emg_load_mendeley:badShape', ...
          'Expected 4 columns, got %d.  Is this the raw file?', size(x,2));
end
N = size(x,1);
t = (0:N-1)'/fs;
fprintf('\n  %s\n  %d samples, %.1f s at %g Hz\n', csvfile, N, N/fs, fs);
if abs(N/fs - 640) > 5
    warning('emg_load_mendeley:duration', ...
        'Duration %.1f s differs from the published 640 s.', N/fs);
end

% ---- quality ---------------------------------------------------------------
Q = emg_quality(x, fs, 'channels', {'col1','col2','col3','col4'}, ...
                'mains', 50, 'mode', 'archival');

% ---- filter and envelope ---------------------------------------------------
Xf  = emg_filter(x, fs, 'band',[20 450], 'mains',50, 'harmonics',3);
w   = round(0.25*fs);
env = sqrt(movmean(Xf.^2, w, 1));
envn = env ./ prctile_local(env, 99);

% ---- inter-channel envelope correlation  (a real diagnostic) ---------------
C = corrcoef(env);
offdiag = C(~eye(4));
fprintf('\n  ENVELOPE CORRELATION between columns\n');
disp(round(C,3));
if max(offdiag) > 0.85
    fprintf(['  >> WARNING: columns correlate at %.2f.  Four DIFFERENT forearm\n' ...
             '     muscles performing flexion vs extension should not track each\n' ...
             '     other this closely.  Likely crosstalk or a shared reference.\n' ...
             '     A 2-channel subset taken from correlated channels carries less\n' ...
             '     independent information than it appears to.  Report this.\n'], ...
             max(offdiag));
end

% ---- segmentation ----------------------------------------------------------
% TWO STAGES.
%   Stage 1: find the four long inter-cycle rests.  These are unambiguous and
%            they anchor the five cycles.
%   Stage 2: INSIDE each cycle, find the individual contractions.  Do NOT lay
%            down a fixed 10 s grid - the earlier version did, and it drifted
%            out of sync by the end of a cycle, which scrambled the labels and
%            made the channel identification fail.  Real contractions are 6 s
%            of activity separated by 4 s of quiet, which is easy to detect.
%
% NOTE the first gesture in each cycle is 'rest', which produces NO burst.
% So a correctly segmented cycle contains NINE bursts, not ten.
tot = mean(envn,2);
dec = 200; te = tot(1:dec:end); fse = fs/dec;
[bs,as] = butter(2, 0.15/(fse/2), 'low');
ts = filtfilt(bs, as, te);
quiet = ts < prctile_local(ts, 35);
[qs, qe] = runs(quiet);
qdur = (qe-qs)/fse;
islong = qdur > 0.6*o.t_cycle;
fprintf('\n  SEGMENTATION\n');
fprintf('    stage 1: long quiet periods found: %d  (expect %d inter-cycle rests)\n', ...
        sum(islong), o.ncycle-1);

nexp = numel(o.gestures) - 1;              % 9 bursts expected per cycle
if sum(islong) >= o.ncycle-1
    cyc_start = [0; qe(islong)/fse];
    cyc_end   = [qs(islong)/fse; N/fs];
    cyc_start = cyc_start(1:o.ncycle);
    cyc_end   = cyc_end(1:o.ncycle);
    src = 'anchored on detected inter-cycle rests';
else
    step = (N/fs)/o.ncycle;
    cyc_start = (0:o.ncycle-1)'*step;
    cyc_end   = cyc_start + step;
    src = 'UNIFORM FALLBACK - inter-cycle rests not cleanly detected';
end
fprintf('    cycle windows (s): %s   [%s]\n', ...
        strjoin(arrayfun(@(a,b)sprintf('%.0f-%.0f',a,b), cyc_start, cyc_end, 'uni',0), '  '), src);

% Stage 2.  IMPORTANT - we cannot simply enumerate bursts.  Of the ten gestures,
% only about five strongly load these four WRIST muscles: extension, flexion,
% ulnar deviation, radial deviation and grip.  Finger abduction/adduction use the
% intrinsic hand muscles, and supination/pronation use supinator and pronator
% teres - none of which are recorded here.  Those four gestures are therefore
% NEARLY SILENT on ECR/FCR/ECU/FCU, and a burst detector finds 4-6 bursts per
% cycle, not 9.  That is the data behaving correctly, not a detection failure.
%
% So: anchor the fixed 10 s protocol grid on the FIRST burst of each cycle.
% Gesture 1 is 'rest' (silent) and gesture 2 is 'extension' (reliably strong),
% so the first burst in a cycle is extension, and everything else follows from
% the published timing.  One strong anchor per cycle beats a rest-boundary
% estimate, and it beats trying to find gestures that are not there.
seg = []; anchors_ok = 0;
t_slot = o.t_hold + o.t_rest;          % 10 s per gesture
for cy = 1:o.ncycle
    i0 = max(1, round(cyc_start(cy)*fs)+1);
    i1 = min(N, round(cyc_end(cy)*fs));
    if i1 - i0 < 20*fs, continue; end
    e = tot(i0:i1);

    lo = prctile_local(e, 15);
    hi = prctile_local(e, 98);
    th = lo + 0.45*(hi - lo);
    [bsm,asm] = butter(2, 0.8/(fs/2), 'low');
    sm = filtfilt(bsm, asm, e);
    [bst, bsp] = runs(sm > th);
    keep = (bsp - bst)/fs > 2.0;
    bst = bst(keep);  bsp = bsp(keep);

    if isempty(bst)
        fprintf('    cycle %d: no anchor burst found - SKIPPED\n', cy);
        continue;
    end
    % first burst = gesture 2 (extension).  Its HOLD starts at the burst onset.
    t_ext = (i0 + bst(1) - 1)/fs;
    t_g1  = t_ext - t_slot;                    % back up one slot to 'rest'
    anchors_ok = anchors_ok + 1;
    fprintf('    cycle %d: anchor (extension onset) at %6.1f s  -> grid starts %6.1f s\n', ...
            cy, t_ext, t_g1);

    for g = 1:numel(o.gestures)
        t0 = t_g1 + (g-1)*t_slot + 0.5;        % skip the first 0.5 s of onset
        t1 = t0 + o.t_hold - 1.0;              % and the last 0.5 s of release
        j0 = round(t0*fs)+1;  j1 = round(t1*fs);
        if j0 >= 1 && j1 <= N && j1 > j0
            seg(end+1,:) = [j0 j1 g];          %#ok<AGROW>
        end
    end
end
fprintf('    %d of %d cycles anchored\n', anchors_ok, o.ncycle);
fprintf('    %d gesture windows built\n', size(seg,1));

% ---- per-gesture activity profile -----------------------------------------
prof = zeros(numel(o.gestures), 4);
cnt  = zeros(numel(o.gestures), 1);
for k = 1:size(seg,1)
    g = seg(k,3);
    prof(g,:) = prof(g,:) + mean(envn(seg(k,1):seg(k,2), :), 1);
    cnt(g)    = cnt(g) + 1;
end
prof = prof ./ max(cnt,1);
if cnt(1) == 0
    warning('emg_load_mendeley:noRest', ...
        'No rest windows found - activity is not referenced to rest.');
end

% subtract the rest row so we score CHANGE FROM REST, not absolute level
irest = find(strcmpi(o.gestures,'rest'), 1);
act   = prof - prof(irest,:);

fprintf('\n  MEAN ACTIVITY ABOVE REST, gesture x column\n');
fprintf('    %-12s %8s %8s %8s %8s\n','gesture','col1','col2','col3','col4');
for g = 1:numel(o.gestures)
    fprintf('    %-12s %8.3f %8.3f %8.3f %8.3f\n', o.gestures{g}, act(g,:));
end

% ---- ALIGNMENT SELF-CHECK --------------------------------------------------
% If the grid is aligned, EXTENSION and FLEXION must load DIFFERENT channels.
% If the same channel is top for both, the labels are shifted and nothing
% downstream can be trusted.
[~, top_ext] = max(act(strcmpi(o.gestures,'extension'), :));
[~, top_flx] = max(act(strcmpi(o.gestures,'flexion'),   :));
fprintf('\n  ALIGNMENT CHECK\n');
fprintf('    extension loads col%d most;  flexion loads col%d most\n', top_ext, top_flx);
if top_ext == top_flx
    fprintf('    >> FAIL: the same column tops both.  Grid is misaligned.\n');
    aligned = false;
else
    fprintf('    >> PASS: antagonists separate, as they must.\n');
    aligned = true;
end

% ---- CHANNEL IDENTIFICATION ------------------------------------------------
% Expected signature, rows = muscle, cols = [flexion extension radial ulnar]
muscles = {'FCR','FCU','ECR','ECU'};
E = [ 1 0 1 0 ;      % FCR : flexion + radial
      1 0 0 1 ;      % FCU : flexion + ulnar
      0 1 1 0 ;      % ECR : extension + radial
      0 1 0 1 ];     % ECU : extension + ulnar

gi = [find(strcmpi(o.gestures,'flexion')), find(strcmpi(o.gestures,'extension')), ...
      find(strcmpi(o.gestures,'radial_dev')), find(strcmpi(o.gestures,'ulnar_dev'))];
M = act(gi, :)';                      % 4 columns x 4 gestures
M = M - mean(M,2);                    % centre per column
sd = std(M,0,2); sd(sd==0) = 1;
M = M ./ sd;                          % z-score each column's gesture profile
Ez = E - mean(E,2);  Ez = Ez ./ std(Ez,0,2);

P = perms(1:4);
score = zeros(size(P,1),1);
for k = 1:size(P,1)
    score(k) = sum(sum(M .* Ez(P(k,:), :)));   % column j -> muscle P(k,j)
end
[ss, order] = sort(score, 'descend');
best = P(order(1), :);
margin = (ss(1)-ss(2)) / max(abs(ss(1)), eps);

fprintf('\n  CHANNEL IDENTIFICATION  (from physiology, not from a forum post)\n');
fprintf('    global assignment score %.3f, runner-up %.3f, margin %.0f %%\n', ...
        ss(1), ss(2), 100*margin);

% PER-COLUMN DIAGNOSTIC.
% The assignment above is PERMUTATION-CONSTRAINED: each muscle is used exactly
% once, because there is exactly one of each on the arm.  That constraint is
% not optional and must not be discarded.
%
% Separately, we ask of each column: independent of the permutation, which
% muscle does it MOST look like, and by how much?  Where that free choice
% AGREES with the constrained assignment and the margin is healthy, the column
% is trustworthy.  Where it disagrees, or the margin is small, it is not - and
% that is usually crosstalk, not a coding error.
fprintf('\n    %-7s', 'column');
for i = 1:4, fprintf('%8s', muscles{i}); end
fprintf('   assigned  freebest  margin   status\n');
conf = zeros(1,4);  agree = false(1,4);
for j = 1:4
    sc = zeros(1,4);
    for i = 1:4, sc(i) = M(j,:) * Ez(i,:)'; end
    [sorted, ord] = sort(sc, 'descend');
    conf(j)  = (sorted(1)-sorted(2)) / max(abs(sorted(1)), eps);
    agree(j) = strcmp(muscles{ord(1)}, muscles{best(j)});
    if agree(j) && conf(j) >= 0.35
        st = 'TRUST';
    elseif ~agree(j)
        st = 'disagrees';
    else
        st = 'weak';
    end
    fprintf('    col%-4d', j);
    fprintf('%8.2f', sc);
    fprintf('   %-8s  %-8s %5.0f %%   %s\n', ...
            muscles{best(j)}, muscles{ord(1)}, 100*conf(j), st);
    map(j).muscle     = muscles{best(j)};     % permutation-constrained ANSWER
    map(j).free_best  = muscles{ord(1)};      % unconstrained favourite
    map(j).confidence = conf(j);
    map(j).agree      = agree(j);
end

iF = find(strcmp({map.muscle},'FCR'), 1);
iE = find(strcmp({map.muscle},'ECR'), 1);

fprintf('\n    VERDICT\n');
if ~isempty(iF) && ~isempty(iE) && aligned && ...
   conf(iF) >= 0.35 && conf(iE) >= 0.35 && agree(iF) && agree(iE)
    fprintf('      The Mode A pair is well determined:\n');
    fprintf('        FCR = col %d  (%.0f %% confidence)\n', iF, 100*conf(iF));
    fprintf('        ECR = col %d  (%.0f %% confidence)\n', iE, 100*conf(iE));
    fprintf('      PROCEED with the 2-channel montage.\n');
    if any(conf < 0.35)
        w = find(conf < 0.35);
        fprintf(['      Note: col%s are weakly determined and should NOT be\n' ...
                 '      used as independent muscles.  Say so in the write-up.\n'], ...
                 strjoin(arrayfun(@(v)sprintf(' %d',v), w, 'uni',0), ' and'));
    end
    map_ok = true;
else
    fprintf('      Mapping not trustworthy.  Use the dataset''s own\n');
    fprintf('      EMG_gesture_segmentation.m and re-run.\n');
    map_ok = false;
end

iF = find(strcmp({map.muscle},'FCR'));
iE = find(strcmp({map.muscle},'ECR'));
fprintf('\n  MODE A MONTAGE (D-13):  FCR = col %d,  ECR = col %d\n', iF, iE);
fprintf('    two_ch = D.Xf(:, D.mode_a);   %% agonist ECR first for extension\n');

% ---- output ----------------------------------------------------------------
D.x = x; D.Xf = Xf; D.env = env; D.envn = envn; D.t = t; D.fs = fs;
D.seg = seg; D.gestures = o.gestures; D.profile = prof; D.act = act;
D.map = map; D.map_margin = margin; D.aligned = aligned;
D.map_ok = map_ok; D.conf = conf; D.agree = agree;
D.n_anchored = anchors_ok; D.Q = Q; D.corr = C;
D.mode_a = [iE iF];   % [agonist antagonist] for Mode A extension
D.source = csvfile; D.segmentation_source = src;

if o.plot
    figure('Name','Mendeley EMG - envelopes and gesture windows','Color','w');
    for c = 1:4
        subplot(4,1,c); plot(t, envn(:,c), 'LineWidth',0.8); hold on; grid on; box on;
        yl = ylim;
        for k = 1:size(seg,1)
            patch(t([seg(k,1) seg(k,2) seg(k,2) seg(k,1)]), [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.85 0.9 0.95], 'EdgeColor','none', 'FaceAlpha',0.5);
        end
        plot(t, envn(:,c), 'LineWidth',0.8);
        ylabel(sprintf('col%d (%s)', c, map(c).muscle));
        if c==1, title('Normalised RMS envelope, shaded = gesture hold windows'); end
        if c==4, xlabel('time [s]'); end
    end
end
fprintf('\n');
end

% =============================================================================
function s = tern(tf, a, b)
if tf, s = a; else, s = b; end
end

% -----------------------------------------------------------------------------
function [s,e] = runs(mask)
d = diff([false; mask(:); false]);
s = find(d==1); e = find(d==-1)-1;
end

function v = prctile_local(X, p)
X = sort(X,1);  n = size(X,1);
v = X(max(1,min(n,ceil(p/100*n))), :);
end

function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_load_mendeley:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
