function [tau, info] = emg_baseline(emg_ago, emg_ant, fs, varargin)
%EMG_BASELINE  Proportional myoelectric intent estimator - the BASELINE (D-21).
%
%   [tau, info] = EMG_BASELINE(emg_ago, emg_ant, fs)
%   [tau, info] = EMG_BASELINE(..., 'mvc_ago', a, 'mvc_ant', b, 'gain', g)
%
%   This is deliberately the SIMPLEST thing that can work: rectify, envelope,
%   normalise, subtract antagonist, deadband, gain.  Classical proportional
%   myoelectric control - 1970s technology, about forty lines.
%
%   WHY IT EXISTS.  The assistive mode needs to know the patient is TRYING, so
%   that assistance is contingent on effort rather than constant (otherwise the
%   patient slacks - step 1 of the design chain, D-02).  EMG is load-bearing for
%   that.  Machine learning is NOT: this function does the job without it.
%
%   Its role is to be the CONTROL CONDITION.  The autoencoder + classifier
%   pipeline is evaluated against this, not assumed better than it.  Then:
%       ML beats baseline  -> a demonstrated result, with a proper control
%       ML does not        -> a negative result consistent with Wang 2026,
%                             honestly reported, and publishable
%   Either outcome is a thesis chapter.  Neither outcome breaks the device.
%   See DECISION_LOG D-21 and FAILURE_LOG F-05.
%
%   INPUTS
%     emg_ago, emg_ant   raw EMG, agonist and antagonist, same length, volts
%                        Mode A (FE):  ago = ECR (extension), ant = FCR
%                        Mode B (RUD): ago = FCR (radial),    ant = FCU
%     fs                 sample rate of the EMG record, Hz (MP36 rate)
%
%   OUTPUT
%     tau                signed intent torque, N.m, at the EMG rate
%                        positive = agonist direction
%     info               struct with envelopes, activation, and the settings used
%
%   To produce the 50 Hz five-field packet for the replay harness, resample tau
%   to 50 Hz and pair it with a confidence and a valid flag - see the interface
%   contract in docs/CONTROL_HANDBOOK.html Part 6.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev D 2026-09-20

% ---- options ---------------------------------------------------------------
opt.bp        = [20 450];   % Hz    Same band as emg_filter.m; see that file's
                            %       header for the corrected citation (De Luca
                            %       et al. 2010, not "SENIAM" - DECISION_LOG D-25)
opt.notch     = 50;         % Hz    mains; set [] to disable
opt.env_lp    = 3.0;        % Hz    envelope cutoff.  Voluntary wrist movement is
                            %       under ~1 Hz (Interview 1), so 3 Hz passes the
                            %       intent and rejects the rest.
opt.mvc_ago   = [];         % V     normalisation reference, agonist
opt.mvc_ant   = [];         % V     normalisation reference, antagonist
opt.deadband  = 0.08;       % -     fraction of MVC below which intent is zero.
                            %       Rejects resting activity and crosstalk.
opt.gain      = 0.60;       % N.m   torque at full normalised activation.
                            %       0.60 of the 0.98 N.m rated - leaves headroom.
opt.tau_max   = 0.80;       % N.m   hard clamp, matches wrist_params.tau_max
opt.rest_win  = [];         % [i1 i2] index range of a quiet rest period, used to
                            %       estimate the noise floor when MVC is absent
opt = parseopts(opt, varargin);

emg_ago = emg_ago(:);  emg_ant = emg_ant(:);
if numel(emg_ago) ~= numel(emg_ant)
    error('emg_baseline:lengthMismatch', 'Channels must be the same length.');
end

% ---- 1. band-pass ----------------------------------------------------------
[b, a] = butter(4, opt.bp/(fs/2), 'bandpass');
x_ago = filtfilt(b, a, emg_ago);
x_ant = filtfilt(b, a, emg_ant);

% ---- 2. mains notch --------------------------------------------------------
if ~isempty(opt.notch)
    w0 = opt.notch/(fs/2);
    [bn, an] = butter(2, [w0-0.02 w0+0.02], 'stop');
    x_ago = filtfilt(bn, an, x_ago);
    x_ant = filtfilt(bn, an, x_ant);
end

% ---- 3. rectify and envelope ----------------------------------------------
[be, ae] = butter(2, opt.env_lp/(fs/2), 'low');
env_ago = filtfilt(be, ae, abs(x_ago));
env_ant = filtfilt(be, ae, abs(x_ant));

% ---- 4. normalisation ------------------------------------------------------
% Prefer a measured MVC.  Failing that, use a high percentile of the record,
% which is crude but stable, and SAY SO in the thesis rather than implying MVC.
if isempty(opt.mvc_ago), opt.mvc_ago = pct95(env_ago); end
if isempty(opt.mvc_ant), opt.mvc_ant = pct95(env_ant); end
if opt.mvc_ago <= 0 || opt.mvc_ant <= 0
    error('emg_baseline:badMVC', 'Normalisation reference must be positive.');
end
a_ago = env_ago / opt.mvc_ago;
a_ant = env_ant / opt.mvc_ant;

% ---- 5. co-contraction-aware difference ------------------------------------
% The signed difference, not the agonist alone.  A patient co-contracting hard is
% not trying to move - the difference correctly returns near zero, where using the
% agonist channel alone would command assistance into a braced limb.
act = a_ago - a_ant;

% ---- 6. deadband -----------------------------------------------------------
db = opt.deadband;
act_db = sign(act) .* max(abs(act) - db, 0) / max(1 - db, eps);
act_db = max(min(act_db, 1), -1);

% ---- 7. gain and clamp -----------------------------------------------------
tau = opt.gain * act_db;
tau = max(min(tau, opt.tau_max), -opt.tau_max);

% ---- info ------------------------------------------------------------------
info = struct();
info.env_ago = env_ago;  info.env_ant = env_ant;
info.a_ago   = a_ago;    info.a_ant   = a_ant;
info.act     = act;      info.act_db  = act_db;
info.opt     = opt;      info.fs      = fs;

% Rest-based noise floor, if a quiet window was given.  This is the number that
% tells you whether the deadband is set sensibly for THIS subject.
if ~isempty(opt.rest_win)
    r = opt.rest_win(1):opt.rest_win(2);
    info.rest_act_rms = sqrt(mean(act(r).^2));
    info.deadband_over_rest = db / max(info.rest_act_rms, eps);
end

if nargout == 0 || (nargout >= 2 && isfield(info,'rest_act_rms'))
    fprintf('\n  emg_baseline  -  proportional myoelectric intent\n');
    fprintf('  ------------------------------------------------------\n');
    fprintf('  band %g-%g Hz, envelope %.1f Hz, deadband %.2f, gain %.2f N.m\n', ...
            opt.bp(1), opt.bp(2), opt.env_lp, opt.deadband, opt.gain);
    fprintf('  normalisation: ago %.4g V, ant %.4g V\n', opt.mvc_ago, opt.mvc_ant);
    fprintf('  intent range:  %.3f to %.3f N.m\n', min(tau), max(tau));
    if isfield(info,'rest_act_rms')
        fprintf('  rest activation RMS %.3f -> deadband is %.1fx the floor %s\n', ...
                info.rest_act_rms, info.deadband_over_rest, ...
                tern(info.deadband_over_rest >= 3, 'OK', '<-- TOO LOW'));
    end
    fprintf('  ------------------------------------------------------\n');
    fprintf('  This is the CONTROL CONDITION.  Report ML performance\n');
    fprintf('  against it, never on its own.\n\n');
end
end

% =============================================================================
function opt = parseopts(opt, args)
for k = 1:2:numel(args)
    if ~isfield(opt, args{k})
        error('emg_baseline:badOption', 'Unknown option "%s".', args{k});
    end
    opt.(args{k}) = args{k+1};
end
end

function s = tern(tf, a, b)
if tf, s = a; else, s = b; end
end

function v = pct95(x)
% 95th percentile without the Statistics Toolbox.
x = sort(x(:));
n = numel(x);
if n == 0, v = 0; return; end
k = max(1, min(n, ceil(0.95*n)));
v = x(k);
end
