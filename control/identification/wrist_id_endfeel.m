function R = wrist_id_endfeel(theta_deg, tau_Nm, varargin)
%WRIST_ID_ENDFEEL  TEST 2b - put a number on "end-feel".
%
%   R = WRIST_ID_ENDFEEL(theta_deg, tau_Nm)
%   R = WRIST_ID_ENDFEEL(..., 'terminal_frac', 0.20, 'plot', true)
%
%   NO NEW EXPERIMENT REQUIRED.  This runs on the Test 2 quasi-static ramp data
%   you are already collecting.
%
%   ---------------------------------------------------------------------------
%   WHY THIS EXISTS
%
%   Interview 2, on the Modified Ashworth Scale:
%
%     "It relies heavily on tactile feedback ('end-feel') that the therapist must
%      feel... It can vary slightly between physios - what feels stiff to me might
%      feel less stiff to a male physio with more physical strength."
%
%   That is a practitioner stating the MECHANISM of MAS inter-rater unreliability:
%   the grade depends on the examiner's own strength.  Our device applies a known
%   torque at a known velocity, so examiner strength leaves the measurement.
%
%   And "end-feel" has a direct engineering translation.  It is the SHAPE of the
%   torque-angle curve near end of range: a hard end-feel is dtau/dtheta rising
%   steeply, a soft end-feel is a gradual rise.  This function measures that.
%
%   OUTPUTS
%     .K_mid        stiffness over the middle of the range, N.m/rad
%     .K_end        stiffness over the terminal fraction, N.m/rad
%     .endfeel_idx  K_end / K_mid  - how much the curve hardens toward the end
%                     ~1      linear spring, no distinct end-feel
%                     2 - 4   firm
%                     > 4     hard / abrupt
%     .theta_knee   angle at which dtau/dtheta first exceeds the knee threshold
%     .hyst_area    hysteresis loop area, N.m.rad (if both directions supplied)
%
%   INPUTS
%     theta_deg, tau_Nm   one ramp direction, monotonic in theta, from Test 2.
%                         Pass a single sweep; call twice for the two directions.
%
%   Stiffness is NOT symmetric about neutral - fit each side separately, as
%   wrist_id_static.m already does.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev D 2026-09-20

opt.terminal_frac = 0.20;   % last 20% of the tested range = "terminal"
opt.mid_band      = [0.25 0.65];  % fraction of range used for K_mid
opt.knee_mult     = 2.0;    % knee = where dtau/dtheta exceeds knee_mult * K_mid
opt.smooth_n      = 15;     % samples for the derivative smoother
opt.plot          = true;
opt = parseopts(opt, varargin);

th  = deg2rad(theta_deg(:));
tau = tau_Nm(:);
if numel(th) ~= numel(tau)
    error('wrist_id_endfeel:lengthMismatch','theta and tau must match.');
end
if numel(th) < 50
    error('wrist_id_endfeel:tooFewSamples', ...
          'Need >= 50 samples across the ramp; got %d.', numel(th));
end

% ---- order by angle and define the range ----------------------------------
[th, idx] = sort(th);   tau = tau(idx);
th0 = th(1);  th1 = th(end);
span = th1 - th0;
if span <= 0
    error('wrist_id_endfeel:noRange','theta has no range.');
end
f = (th - th0) / span;          % 0..1 through the sweep

% ---- K_mid : straight-line fit over the middle band -----------------------
m = f >= opt.mid_band(1) & f <= opt.mid_band(2);
if nnz(m) < 10
    error('wrist_id_endfeel:midBandEmpty','Middle band has too few samples.');
end
pm = [th(m) ones(nnz(m),1)] \ tau(m);
K_mid = pm(1);

% ---- refine: K_mid must be fitted BELOW the knee, or the stiffening region
%      contaminates it.  Two passes: provisional K_mid -> find knee -> refit
%      K_mid strictly below the knee with a 10% margin -> re-find the knee.
nsm = max(3, round(opt.smooth_n));
for pass = 1:2
    dtau_p = gradient(movmean(tau, nsm), movmean(th, nsm));
    ik_p = find(dtau_p > opt.knee_mult*K_mid & f > opt.mid_band(1), 1, 'first');
    if isempty(ik_p), break; end
    f_knee = f(ik_p);
    m_new  = f >= opt.mid_band(1) & f <= max(0.9*f_knee, opt.mid_band(1)+0.05);
    if nnz(m_new) < 10, break; end
    pm = [th(m_new) ones(nnz(m_new),1)] \ tau(m_new);
    K_mid = pm(1);
    m = m_new;
end

% ---- K_end : straight-line fit over the terminal fraction -----------------
e = f >= (1 - opt.terminal_frac);
if nnz(e) < 10
    error('wrist_id_endfeel:terminalEmpty','Terminal band has too few samples.');
end
pe = [th(e) ones(nnz(e),1)] \ tau(e);
K_end = pe(1);

endfeel_idx = K_end / K_mid;

% ---- local stiffness curve and the knee -----------------------------------
n = max(3, round(opt.smooth_n));
tau_s = movmean(tau, n);
th_s  = movmean(th,  n);
dtau  = gradient(tau_s, th_s);       % local dtau/dtheta, N.m/rad

thresh = opt.knee_mult * K_mid;
ik = find(dtau > thresh & f > opt.mid_band(2), 1, 'first');
if isempty(ik)
    theta_knee = NaN;
else
    theta_knee = rad2deg(th(ik));
end

% ---- classify --------------------------------------------------------------
if      endfeel_idx < 1.5, cls = 'soft / near-linear';
elseif  endfeel_idx < 2.5, cls = 'firm';
elseif  endfeel_idx < 4.0, cls = 'hard';
else,                      cls = 'abrupt / bony-hard';
end

R = struct('K_mid',K_mid,'K_end',K_end,'endfeel_idx',endfeel_idx, ...
           'theta_knee_deg',theta_knee,'class',cls, ...
           'range_deg',[rad2deg(th0) rad2deg(th1)], ...
           'dtau_dtheta',dtau,'theta_rad',th,'tau_Nm',tau, ...
           'opt',opt,'run_at',datestr(now,31));

% ---- report ----------------------------------------------------------------
fprintf('\n  TEST 2b  -  end-feel index\n');
fprintf('  ---------------------------------------------------------------\n');
fprintf('  range tested   = %6.1f to %6.1f deg\n', R.range_deg(1), R.range_deg(2));
fprintf('  K_mid          = %7.3f N.m/rad   (middle %.0f-%.0f%% of range)\n', ...
        K_mid, 100*opt.mid_band(1), 100*opt.mid_band(2));
fprintf('  K_end          = %7.3f N.m/rad   (terminal %.0f%%)\n', ...
        K_end, 100*opt.terminal_frac);
fprintf('  END-FEEL INDEX = %7.2f           [%s]\n', endfeel_idx, cls);
if ~isnan(theta_knee)
fprintf('  knee at        = %7.1f deg        (dtau/dtheta > %.1f x K_mid)\n', ...
        theta_knee, opt.knee_mult);
else
fprintf('  knee           =    none detected within the tested range\n');
end
fprintf('  ---------------------------------------------------------------\n');
fprintf('  Report this ALONGSIDE the MAS grade, never instead of it - the\n');
fprintf('  comparison between the two is the point.  (D-19)\n\n');

% ---- plot ------------------------------------------------------------------
if opt.plot
    figure('Name','Test 2b - end-feel','Color','w');
    subplot(2,1,1); hold on; grid on; box on;
    plot(rad2deg(th), tau, '-', 'LineWidth',1.3, 'DisplayName','measured');
    plot(rad2deg(th(m)), th(m)*pm(1)+pm(2), '--', 'LineWidth',1.4, ...
         'DisplayName',sprintf('K_{mid} = %.3f', K_mid));
    plot(rad2deg(th(e)), th(e)*pe(1)+pe(2), '--', 'LineWidth',1.4, ...
         'DisplayName',sprintf('K_{end} = %.3f', K_end));
    if ~isnan(theta_knee), xline(theta_knee, ':', 'knee', 'LineWidth',1.2); end
    ylabel('\tau  [N.m]');
    title(sprintf('End-feel index = %.2f  [%s]', endfeel_idx, cls));
    legend('Location','northwest');

    subplot(2,1,2); hold on; grid on; box on;
    plot(rad2deg(th), dtau, '-', 'LineWidth',1.3);
    yline(K_mid, '--', 'K_{mid}');
    yline(thresh, ':', sprintf('%.1f x K_{mid}', opt.knee_mult));
    xlabel('\theta  [deg]'); ylabel('d\tau/d\theta  [N.m/rad]');
end
end

% =============================================================================
function opt = parseopts(opt, args)
for k = 1:2:numel(args)
    if ~isfield(opt, args{k})
        error('wrist_id_endfeel:badOption','Unknown option "%s".', args{k});
    end
    opt.(args{k}) = args{k+1};
end
end
