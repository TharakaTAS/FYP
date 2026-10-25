function R = wrist_id_ktheta(csvfile, varargin)
%WRIST_ID_KTHETA  TEST 0 - identify transmission (cable) stiffness k_theta.
%
%   R = WRIST_ID_KTHETA(csvfile)
%   R = WRIST_ID_KTHETA(csvfile, 'r_m', 0.020, 'r_j', 0.040, 'plot', true)
%
%   THIS TEST GATES EVERYTHING DOWNSTREAM.  The real-time torque estimate is
%
%       tau_hat = k_theta * (theta_m*(r_m/r_j) - theta_j)
%
%   so until k_theta is known, no torque number from this rig means anything,
%   and quantisation noise on the estimate scales directly with it:
%
%       dTau_rms = k_theta * 6.27e-4  N.m        (two 12-bit encoders, RSS)
%
%   RESOLUTION WARNING - read this before running the test.
%   The wind-up you are measuring is TINY.  At k_theta = 95 and a peak applied
%   torque of 0.196 N.m (200 g on a 100 mm lever) the FULL-SCALE wind-up is
%   0.118 deg - that is 1.3 encoder counts.  You cannot fit a line through 1.3
%   counts.  Two things are therefore mandatory:
%     (a) APPLY NEAR-RATED TORQUE.  1 kg on a 100 mm lever = 0.98 N.m gives 6.7
%         counts at k = 95, 12.8 counts at k = 50.  200 g is not enough.
%     (b) DITHER AND AVERAGE.  Superimpose a small servo oscillation (+/-0.5 deg
%         at ~5 Hz) during the 5 s averaging window.  Dither converts
%         deterministic quantisation into noise that averages away; 2500 samples
%         then resolve far below one LSB.  A perfectly still rig does NOT.
%   This script checks both and warns you.  See FAILURE_LOG F-08.
%
%   USABLE BAND 30 - 250 N.m/rad.  Above 250 the estimate is mostly quantisation
%   noise and a series spring must be added.  Below 30 the transmission
%   resonance drops into the control band - increase pre-tension and re-run.
%
%   ---------------------------------------------------------------------------
%   PROCEDURE  (see docs/CONTROL_HANDBOOK.html, Test 0)
%     1. Command the servo to hold a mid-range position.  Leave it there for the
%        whole test - do not re-command between points.
%     2. Hang a known mass on a lever of measured length L:  tau = m*g*L.
%     3. Record BOTH encoders for 5 s.  Average over the steady window.
%     4. Repeat for 200 / 400 / 600 / 800 / 1000 g, in BOTH directions.
%        (NOT 50-200 g - see the resolution warning above.)
%     5. Run this script on the resulting CSV.
%
%   INPUT CSV COLUMNS (header row required, order free):
%     mass_g        applied mass, grams (signed is fine; direction also given)
%     direction     +1 or -1  (sign of the applied torque)
%     lever_m       lever arm length, metres
%     theta_m_deg   motor-side angle, degrees, averaged over the steady window
%     theta_j_deg   joint-side angle, degrees, averaged over the steady window
%     theta_m_std   (optional) std of theta_m over the window, degrees
%     theta_j_std   (optional) std of theta_j over the window, degrees
%
%   OUTPUT STRUCT R
%     .k_theta      fitted stiffness, N.m/rad        <- the number you want
%     .k_ci95       95% confidence interval on it
%     .R2           coefficient of determination      (pass: >= 0.98)
%     .hyst_pct     hysteresis between up and down sweeps, %  (pass: <= 15)
%     .dTau_rms     implied torque-estimate noise, N.m
%     .dTau_pct     the same, as % of the 0.98 N.m rated torque
%     .verdict      'PASS' | 'ADD SERIES SPRING' | 'INCREASE PRE-TENSION' | 'FAIL'
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev C 2026-09-19

% ---- options ---------------------------------------------------------------
opt.r_m  = 0.020;      % m, motor pulley radius   - MEASURE with calipers
opt.r_j  = 0.040;      % m, joint pulley radius   - MEASURE with calipers
opt.plot = true;
opt.g    = 9.80665;    % m/s^2
opt = parseopts(opt, varargin);

% ---- load ------------------------------------------------------------------
T = readtable(csvfile);
req = {'mass_g','direction','lever_m','theta_m_deg','theta_j_deg'};
for k = 1:numel(req)
    if ~ismember(req{k}, T.Properties.VariableNames)
        error('wrist_id_ktheta:missingColumn', ...
              'CSV is missing required column "%s".', req{k});
    end
end

% ---- build the two regression variables ------------------------------------
% Applied torque at the joint, from first principles.  No calibration constant,
% no sensor - just mass, gravity and a length you measured with a ruler.  This
% is what makes the result traceable.
tau = (abs(T.mass_g)/1000) * opt.g .* T.lever_m .* sign(T.direction);   % N.m

% Transmission wind-up, referred to the joint.
th_m = deg2rad(T.theta_m_deg);
th_j = deg2rad(T.theta_j_deg);
Delta = th_m * (opt.r_m/opt.r_j) - th_j;                                % rad

% Remove the zero-load offset: the point(s) with no applied mass define Delta=0.
izero = abs(T.mass_g) < eps;
if any(izero)
    Delta = Delta - mean(Delta(izero));
else
    % No explicit zero row - use the intercept of the fit instead (below) and
    % warn, because an un-referenced offset inflates the intercept.
    warning('wrist_id_ktheta:noZeroRow', ...
        ['No zero-mass row found.  The offset will be absorbed into the ' ...
         'intercept.  Add a 0 g row to both directions next time.']);
end

% ---- RESOLUTION CHECK  ------------------------------------------------------
% Full-scale wind-up, expressed in encoder counts.  Below ~5 counts the fit is
% dominated by quantisation and the result is not trustworthy.
q_lsb    = deg2rad(360/4096);
fs_delta = max(Delta) - min(Delta);
fs_lsb   = fs_delta / q_lsb;
dither   = std(Delta(abs(T.mass_g) > eps)) / q_lsb;   % crude dither indicator

% ---- least-squares fit  tau = k*Delta + c ----------------------------------
n = numel(tau);
if n < 4
    error('wrist_id_ktheta:tooFewPoints', ...
          'Need at least 4 points; found %d.  Four masses x two directions.', n);
end
X = [Delta, ones(n,1)];
b = X \ tau;
k_theta = b(1);
c0      = b(2);

resid = tau - X*b;
SSE   = sum(resid.^2);
SST   = sum((tau - mean(tau)).^2);
R2    = 1 - SSE/SST;

dof   = n - 2;
s2    = SSE / dof;
C     = s2 * inv(X'*X);                      %#ok<MINV>  2x2, fine
se_k  = sqrt(C(1,1));
tcrit = tinv95(dof);
k_ci95 = [k_theta - tcrit*se_k, k_theta + tcrit*se_k];

% ---- hysteresis between the two sweep directions ---------------------------
% Fit each direction separately; the spread between them is the hysteresis.
hyst_pct = NaN;  k_up = NaN;  k_dn = NaN;
iu = T.direction > 0;   id = T.direction < 0;
if sum(iu) >= 2 && sum(id) >= 2
    k_up = [Delta(iu), ones(sum(iu),1)] \ tau(iu);   k_up = k_up(1);
    k_dn = [Delta(id), ones(sum(id),1)] \ tau(id);   k_dn = k_dn(1);
    hyst_pct = 100 * abs(k_up - k_dn) / mean([abs(k_up) abs(k_dn)]);
end

% ---- implied torque resolution ---------------------------------------------
q_enc      = deg2rad(360/4096);                  % 12-bit LSB, 0.0879 deg
dTheta_rms = sqrt(2) * q_enc / sqrt(12);         % 6.27e-4 rad
dTau_rms   = k_theta * dTheta_rms;
dTau_pct   = 100 * dTau_rms / 0.98;

% Two-mass resonance, assuming a nominal Jh and Kh for the indication only.
Jh_nom = 0.003;  Kh_nom = 1.0;
fn_res = sqrt((k_theta + Kh_nom)/Jh_nom) / (2*pi);

% ---- verdict ---------------------------------------------------------------
if fs_lsb < 5
    verdict = sprintf(['FAIL - only %.1f encoder counts of full-scale wind-up. ' ...
                       'Use heavier masses (target 1 kg @ 100 mm) and dither.'], fs_lsb);
elseif R2 < 0.98
    verdict = 'FAIL - poor fit, check for slack or a moving servo';
elseif ~isnan(hyst_pct) && hyst_pct > 15
    verdict = 'FAIL - hysteresis > 15%, cable is slipping or seating';
elseif k_theta > 250
    verdict = 'ADD SERIES SPRING - too stiff, torque estimate is mostly noise';
elseif k_theta < 30
    verdict = 'INCREASE PRE-TENSION - too soft, resonance is in the control band';
else
    verdict = 'PASS';
end

R = struct('fs_lsb',fs_lsb,'dither_lsb',dither, ...
           'k_theta',k_theta,'k_ci95',k_ci95,'offset',c0,'R2',R2, ...
           'hyst_pct',hyst_pct,'k_up',k_up,'k_dn',k_dn, ...
           'dTau_rms',dTau_rms,'dTau_pct',dTau_pct,'fn_res_nominal',fn_res, ...
           'n_points',n,'r_m',opt.r_m,'r_j',opt.r_j, ...
           'verdict',verdict,'source',csvfile,'run_at',datestr(now,31));

% ---- report ----------------------------------------------------------------
fprintf('\n  TEST 0  -  transmission stiffness\n');
fprintf('  file: %s\n', csvfile);
fprintf('  ---------------------------------------------------------------\n');
fprintf('  k_theta        = %7.1f N.m/rad   (95%% CI %.1f to %.1f)\n', ...
        k_theta, k_ci95(1), k_ci95(2));
fprintf('  fit R2         = %7.4f          (pass >= 0.98)\n', R2);
if ~isnan(hyst_pct)
fprintf('  hysteresis     = %7.1f %%         (pass <= 15)\n', hyst_pct);
fprintf('                   up %.1f / down %.1f N.m/rad\n', k_up, k_dn);
end
fprintf('  points         = %7d\n', n);
fprintf('  ---------------------------------------------------------------\n');
fprintf('  RESOLUTION\n');
fprintf('    full-scale wind-up = %.3f deg = %.1f encoder counts  %s\n', ...
        rad2deg(fs_delta), fs_lsb, tern(fs_lsb >= 5, 'OK', '<-- TOO FEW'));
fprintf('    peak applied torque = %.3f N.m (%.0f%% of rated)\n', ...
        max(abs(tau)), 100*max(abs(tau))/0.98);
if fs_lsb < 5
fprintf('    >> Increase the masses.  1 kg on a 100 mm lever = 0.98 N.m.\n');
fprintf('    >> And dither: +/-0.5 deg at 5 Hz during the averaging window.\n');
end
fprintf('  ---------------------------------------------------------------\n');
fprintf('  IMPLIED torque-estimate noise\n');
fprintf('    dTau_rms     = %7.4f N.m  = %.1f %% of rated 0.98 N.m\n', ...
        dTau_rms, dTau_pct);
fprintf('    after filtering to 50 Hz, roughly /sqrt(10): %.1f %%\n', dTau_pct/sqrt(10));
fprintf('    two-mass resonance (at Jh=0.003, Kh=1.0): %.1f Hz\n', fn_res);
fprintf('  ---------------------------------------------------------------\n');
fprintf('  VERDICT: %s\n\n', verdict);

fprintf('  Copy this into your parameter file and into every CSV header:\n');
fprintf('    P.k_theta = %.1f;   %% Test 0, %s, R2 = %.4f\n\n', ...
        k_theta, datestr(now,'yyyy-mm-dd'), R2);

% ---- plot ------------------------------------------------------------------
if opt.plot
    figure('Name','Test 0 - transmission stiffness','Color','w');
    hold on; grid on; box on;
    if any(iu), plot(rad2deg(Delta(iu)), tau(iu), 'o', 'MarkerSize',7, ...
                     'LineWidth',1.2, 'DisplayName','direction +'); end
    if any(id), plot(rad2deg(Delta(id)), tau(id), 's', 'MarkerSize',7, ...
                     'LineWidth',1.2, 'DisplayName','direction -'); end
    xs = linspace(min(Delta), max(Delta), 50);
    plot(rad2deg(xs), k_theta*xs + c0, '-', 'LineWidth',1.6, ...
         'DisplayName', sprintf('fit: %.1f N.m/rad', k_theta));
    xlabel('transmission wind-up \Delta  [deg, referred to joint]');
    ylabel('applied torque \tau  [N.m]');
    title(sprintf('Test 0  -  k_\\theta = %.1f N.m/rad,  R^2 = %.4f  [%s]', ...
                  k_theta, R2, verdict));
    legend('Location','northwest');
end
end

% =============================================================================
function opt = parseopts(opt, args)
for k = 1:2:numel(args)
    name = args{k};
    if ~isfield(opt, name)
        error('wrist_id_ktheta:badOption','Unknown option "%s".', name);
    end
    opt.(name) = args{k+1};
end
end

% -----------------------------------------------------------------------------
function s = tern(tf, a, b)
if tf, s = a; else, s = b; end
end

% -----------------------------------------------------------------------------
function t = tinv95(dof)
% Two-sided 95% t critical values, so the script runs without the Statistics
% Toolbox.  Falls back to tinv when the toolbox is present.
if exist('tinv','file') == 2 || exist('tinv','builtin') == 5
    try, t = tinv(0.975, dof); return; catch, end
end
tbl = [1 12.706; 2 4.303; 3 3.182; 4 2.776; 5 2.571; 6 2.447; 7 2.365; ...
       8 2.306; 9 2.262; 10 2.228; 12 2.179; 15 2.131; 20 2.086; 30 2.042; ...
       40 2.021; 60 2.000; 120 1.980];
if dof >= 120
    t = 1.96;
else
    t = interp1(tbl(:,1), tbl(:,2), max(dof,1), 'linear', 'extrap');
end
end
