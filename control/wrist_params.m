function P = wrist_params(condition)
%WRIST_PARAMS  Parameter set for the MC 4206 wrist exoskeleton control design.
%
%   P = WRIST_PARAMS()             returns the 'normal' condition
%   P = WRIST_PARAMS(condition)    condition = 'relaxed' | 'normal' | 'spastic' | 'vspastic'
%
%   Every field carries a STATUS comment:
%       [KNOWN]   established from datasheet or literature
%       [CHOSEN]  our design decision - must be justified in the thesis
%       [MEASURE] placeholder; replace with the identified value
%
%   Rev D - 2026-09-20.  Adds the clinical-requirement block (D-19), the protocol
%   block revised after the physiotherapist interviews (D-17, D-18), and an
%   adjustable device neutral (D-15).  See docs/CONTROL_HANDBOOK.html Part 3 for
%   the full parameter register, docs/DECISION_LOG.md for why each value is what
%   it is, and docs/PHYSIO_FIELD_NOTES.md for what the clinic actually does.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka

if nargin < 1 || isempty(condition)
    condition = 'normal';
end
P.condition = lower(condition);

% ----------------------------------------------------------------------------
% 1.  HUMAN LIMB  -  the plant.  Eq. (2.1) in the handbook.
%     J_h*thdd + B_h*thd + K_h*(th - th0) + tau_spas = tau_dev + tau_grav + tau_vol
% ----------------------------------------------------------------------------
% J_h is dominated by the hand plus the device hand-plate and barely varies with
% impairment, so it is held constant across conditions.  B_h and K_h do vary.
% Literature bracket for K_h: Falzarano et al. 2021, LIVE subjects, wrist FE,
% near neutral: 0.43 - 1.89 N.m/rad.  See FAILURE_LOG F-04 before using any
% larger number you find in the cadaver literature.

P.Jh = 0.003;                       % kg.m^2      [MEASURE] Test 1 / 1b

switch P.condition
    case 'relaxed'
        P.Bh = 0.015;               % N.m.s/rad   [MEASURE] Test 1 / 3
        P.Kh = 0.45;                % N.m/rad     [MEASURE] Test 2 / 3
    case 'normal'
        P.Bh = 0.030;
        P.Kh = 1.00;
    case 'spastic'
        P.Bh = 0.050;
        P.Kh = 1.60;
    case 'vspastic'
        P.Bh = 0.070;
        P.Kh = 2.40;                % deliberately beyond the Falzarano upper
                                    % bound, as a stress case for the margins
    otherwise
        error('wrist_params:badCondition', ...
              'condition must be relaxed | normal | spastic | vspastic');
end

% theta0 is NOT zero in this population.  Interview 2: "In a neutral position the
% hand rests slightly bent... it stays slightly curved downward."  In wrist drop the
% hand rests FLEXED, so the device neutral must be set to the patient's own resting
% posture before strapping and the offset recorded (D-15).  A fixed anatomical-zero
% neutral collapses the usable extension range - the side that matters clinically.
P.theta0 = 0;                       % rad         [MEASURE] Test 2, per subject

% ----------------------------------------------------------------------------
% 2.  SPASTICITY  -  Eq. (2.2).  Lance (1980): velocity-dependent stretch reflex.
%     tau_spas = bs * max(|thd| - thd_th, 0) * sign(thd)
% ----------------------------------------------------------------------------
switch P.condition
    case 'relaxed',  P.bs = 0.00;  P.thd_th = Inf;   % no reflex at all
    case 'normal',   P.bs = 0.00;  P.thd_th = Inf;   % healthy control -> bs = 0
    case 'spastic',  P.bs = 0.12;  P.thd_th = 1.20;  % [MEASURE] Test 4
    case 'vspastic', P.bs = 0.28;  P.thd_th = 0.60;  % [MEASURE] Test 4
end
% NOTE: on a healthy subject Test 4 must return bs ~ 0.  That is the negative
% control and it has to appear in the thesis.

% ----------------------------------------------------------------------------
% 3.  TRANSMISSION  -  CORRECTED at Rev D (D-23).  D-14's "one servo per DOF,
%     cable loop" was never uniform across both axes.  RUD has always been
%     cable-loop driven.  FE has always been direct motor-driven (rigid link,
%     no cable) - D-14 misstated this as uniform.  A direct link has no
%     windup to measure, so FE gets its own added series spring (PHYSICAL_BOM
%     Sec 7.3) to restore the same free two-encoder torque estimate RUD
%     already has:  tau_hat = k_theta * Delta.  The two axes have TWO
%     independent transmission stiffnesses now, not one.
% ----------------------------------------------------------------------------
P.k_theta_rud = 80;                 % N.m/rad     [MEASURE] Test 0, RUD (cable)  <-- GATES RUD
P.k_theta_fe  = NaN;                % N.m/rad     [MEASURE] Test 0, FE (added spring, D-23)
                                    %   Spring constant still being sized by the mechanical
                                    %   team as of D-23 (2026-09-22).  NaN on purpose - do not
                                    %   silently default this to a guessed number.  Once fitted,
                                    %   confirm it lands in the 40-60 N.m/rad target band
                                    %   (PHYSICAL_BOM Sec 7.2) before trusting FE's torque loop.
P.r_m       = 0.020;                % m           [MEASURE] calipers, RUD motor pulley
P.r_j       = 0.040;                % m           [MEASURE] calipers, RUD joint pulley
P.n_ratio   = P.r_m / P.r_j;        % -           derived, RUD only

% Torque-estimate resolution, Eq. (2.10).  Both encoders are 12-bit over 360 deg.
% Quantisation noise is a property of the ENCODERS, so it's the same for both axes;
% the resulting torque noise (dTau_rms) still depends on each axis's own k_theta.
P.q_enc      = deg2rad(360/4096);              % rad, one LSB = 0.0879 deg
P.dTheta_rms = sqrt(2) * P.q_enc / sqrt(12);   % rad, RSS of two quantisers
P.dTau_rms_rud = P.k_theta_rud * P.dTheta_rms; % N.m RMS on RUD's torque estimate
P.dTau_rms_fe  = P.k_theta_fe  * P.dTheta_rms; % N.m RMS on FE's torque estimate (NaN until spring is sized)

% ----------------------------------------------------------------------------
% 4.  ACTUATOR  -  Feetech STS3215, 12 V, 30 kgf.cm.  Eq. (2.11).
% ----------------------------------------------------------------------------
P.fs_servo  = 8.0;                  % Hz          [MEASURE] Test 8 (assumed)
P.zs_servo  = 0.8;                  % -           [MEASURE] Test 8 (assumed)
P.Td        = 0.005;                % s           [MEASURE] Test 7 (assumed)
P.tau_rated = 0.98;                 % N.m         [KNOWN] datasheet
P.tau_stall = 2.94;                 % N.m         [KNOWN] datasheet
P.tau_max   = 0.80;                 % N.m         [CHOSEN] working limit, 82% of rated
P.rate_max  = 0.60;                 % rad/s       [CHOSEN] safety rate limit
P.theta_lim = deg2rad(45);          % rad         [CHOSEN] = mechanical hard stop
P.theta_lim_rom = deg2rad(60);      % rad         [CHOSEN] ROM-work stop setting.
                                    %   Normal wrist ROM is ~80 deg flexion / ~70 deg
                                    %   extension, so +/-45 CANNOT deliver terminal-
                                    %   range stretching.  Stops must be adjustable;
                                    %   45 for measurement and assistance, 60 under
                                    %   supervision for ROM work.  (D-19)
P.temp_stop = 60;                   % degC        [CHOSEN] abort threshold (D-14)

% ----------------------------------------------------------------------------
% 5.  TARGET IMPEDANCE  -  what we render.  Eqs. (2.4) - (2.7).  All [CHOSEN].
% ----------------------------------------------------------------------------
P.Md   = 0.010;                     % kg.m^2      chosen so Md >= 3*Jh
P.fn_d = 2.0;                       % Hz          above the voluntary band (<1 Hz)
P.zeta_d = 0.8;                     % -           settles without overshooting
P.Kd = P.Md * (2*pi*P.fn_d)^2;                   % = 1.5791 N.m/rad
P.Bd = 2 * P.zeta_d * sqrt(P.Kd * P.Md);         % = 0.2011 N.m.s/rad

% Assist-as-Needed levels.  Kd is the assistance knob (D-04).
P.AAN = [0.4 0.8 1.6 3.2 6.4];      % N.m/rad
P.AAN_default = 3;                  % index -> 1.6, closest to the nominal Kd

% ----------------------------------------------------------------------------
% 6.  TIMING  -  multirate architecture, handbook Part 6.2.
% ----------------------------------------------------------------------------
P.Ts_ctrl   = 0.002;                % s   500 Hz  admittance + sync-write + log
P.Ts_torque = 0.005;                % s   200 Hz  servo readback -> torque estimate
P.Ts_intent = 0.020;                % s    50 Hz  ML intent, zero-order held
P.f_lp_tau  = 9.0;                  % Hz  [CHOSEN] roll-off on tau_ext, keeps the
                                    %     transmission resonance out of the loop

% ----------------------------------------------------------------------------
% 6b. CLINICAL REQUIREMENTS  -  derived from practice, not from the datasheet.
%     Interview 2: early-phase goals progress ~1 cm of hand lift per week.
%     At a wrist-axis-to-landmark distance of 8.5 cm (3rd MCP) to 18 cm (fingertip)
%     that is a weekly clinical increment of roughly 3 to 7 degrees.  A factor of
%     three below the smallest of those gives the requirement below.  (D-19)
% ----------------------------------------------------------------------------
P.MDC_target_deg  = 2.0;            % deg  [CHOSEN] minimum detectable change, Test 9
P.MDC_stretch_deg = 1.0;            % deg  stretch goal
P.weekly_change_deg = [3 7];        % deg  the clinical increment we must resolve

% ----------------------------------------------------------------------------
% 6c. PROTOCOL  -  constrained by clinical practice, not by our convenience.
%     Interview 2: "we flex and extend at most twice; we do not exceed that."
%     Sweeping velocity in order across many passes confounds repetition with
%     velocity SYSTEMATICALLY - the bias looks like the effect.  (D-17, F-06)
% ----------------------------------------------------------------------------
P.T4_velocities   = [10 30 60 120 240];   % deg/s   graded stretch velocities
P.T4_passes       = 2;                    % per velocity  (was 3)
P.T4_randomise    = true;                 % velocity order, per subject
P.T4_rest_s       = 60;                   % minimum, logged
P.T4_patient_max  = 2;                    % TOTAL passes in a patient protocol
P.T4b_passes      = 10;                   % Test 4b conditioning decay, healthy only

% Posture is a CONTROLLED variable.  FDS/FDP/EDC cross both wrist and fingers, so
% finger posture changes measured wrist stiffness.  (D-18, F-07)
P.finger_posture      = 'constrained_grip';   % standardised, mechanically held
P.forearm_rotation_deg = 0;                   % fixed and recorded

% ----------------------------------------------------------------------------
% 7.  DERIVED QUANTITIES AND DESIGN-RULE CHECKS  -  now per axis (D-23).
%     RUD's resonance uses its cable k_theta; FE's uses its (still-unsized)
%     spring k_theta and will read NaN until that spring constant is fixed.
% ----------------------------------------------------------------------------
P.wn_res_rud = sqrt((P.k_theta_rud + P.Kh) / P.Jh);   % rad/s, Eq. (2.13), RUD
P.fn_res_rud = P.wn_res_rud / (2*pi);                 % Hz
P.wn_res_fe  = sqrt((P.k_theta_fe  + P.Kh) / P.Jh);   % rad/s, Eq. (2.13), FE
P.fn_res_fe  = P.wn_res_fe  / (2*pi);                 % Hz

P.check.Md_over_Jh    = P.Md / P.Jh;              % want >= 3
P.check.Kd_over_Kh    = P.Kd / P.Kh;              % want >= 0.6
P.check.dTau_percent_rud = 100 * P.dTau_rms_rud / P.tau_rated;
P.check.dTau_percent_fe  = 100 * P.dTau_rms_fe  / P.tau_rated;
P.check.fres_over_fs_rud = P.fn_res_rud / P.fs_servo;  % want >= 2
P.check.fres_over_fs_fe  = P.fn_res_fe  / P.fs_servo;  % want >= 2

% ----------------------------------------------------------------------------
% 8.  REPORT  -  printed whenever the function is called without an output.
% ----------------------------------------------------------------------------
if nargout == 0
    fprintf('\n  wrist_params  -  condition: %s\n', P.condition);
    fprintf('  ---------------------------------------------------------\n');
    fprintf('  PLANT       Jh %.4f kg.m^2   Bh %.3f   Kh %.2f N.m/rad\n', ...
            P.Jh, P.Bh, P.Kh);
    fprintf('  SPASTICITY  bs %.3f          threshold %.2f rad/s\n', ...
            P.bs, P.thd_th);
    fprintf('  TRANSMSSN   RUD (cable) k_theta %.0f N.m/rad   ratio r_m/r_j = %.3f\n', ...
            P.k_theta_rud, P.n_ratio);
    fprintf('              FE (spring, D-23) k_theta %s   %s\n', ...
            tern_str(isnan(P.k_theta_fe), '[not yet sized]', sprintf('%.0f N.m/rad', P.k_theta_fe)), ...
            tern_str(isnan(P.k_theta_fe), '<- mechanical team still calculating', ''));
    fprintf('  TARGET      Md %.4f  Bd %.4f  Kd %.4f  (fn %.1f Hz, z %.2f)\n', ...
            P.Md, P.Bd, P.Kd, P.fn_d, P.zeta_d);
    fprintf('  ---------------------------------------------------------\n');
    fprintf('  CHECK  Md/Jh          = %5.2f   (want >= 3.0)    %s\n', ...
            P.check.Md_over_Jh, pass(P.check.Md_over_Jh >= 3));
    fprintf('  CHECK  Kd/Kh          = %5.2f   (want >= 0.6)    %s\n', ...
            P.check.Kd_over_Kh, pass(P.check.Kd_over_Kh >= 0.6));
    fprintf('  CHECK  RUD torq noise = %5.2f %% of rated       %s\n', ...
            P.check.dTau_percent_rud, pass(P.check.dTau_percent_rud <= 10));
    fprintf('  CHECK  RUD f_res      = %5.1f Hz  (%.1fx servo) %s\n', ...
            P.fn_res_rud, P.check.fres_over_fs_rud, pass(P.check.fres_over_fs_rud >= 2));
    fprintf('  CHECK  FE  torq noise = %s\n', ...
            tern_str(isnan(P.k_theta_fe), '(FE spring not yet sized - cannot check)', ...
                     sprintf('%5.2f %% of rated       %s', P.check.dTau_percent_fe, pass(P.check.dTau_percent_fe <= 10))));
    fprintf('  CHECK  FE  f_res      = %s\n', ...
            tern_str(isnan(P.k_theta_fe), '(FE spring not yet sized - cannot check)', ...
                     sprintf('%5.1f Hz  (%.1fx servo) %s', P.fn_res_fe, P.check.fres_over_fs_fe, pass(P.check.fres_over_fs_fe >= 2))));
    fprintf('  ---------------------------------------------------------\n');
    fprintf('  CLINICAL  MDC target   <= %.1f deg   (Test 9 pass criterion)\n', ...
            P.MDC_target_deg);
    fprintf('            weekly change  %.0f-%.0f deg to resolve (Interview 2)\n', ...
            P.weekly_change_deg(1), P.weekly_change_deg(2));
    fprintf('  PROTOCOL  Test 4: %d passes x %d velocities, randomised, %ds rest\n', ...
            P.T4_passes, numel(P.T4_velocities), P.T4_rest_s);
    fprintf('            patient protocol: %d passes TOTAL\n', P.T4_patient_max);
    fprintf('  POSTURE   fingers %s, forearm %.0f deg - both RECORDED\n', ...
            P.finger_posture, P.forearm_rotation_deg);
    fprintf('  ---------------------------------------------------------\n');
    fprintf('  Fields still marked [MEASURE] are placeholders.  Replace them\n');
    fprintf('  with identified values before quoting any margin in the thesis.\n\n');
end
end

% -----------------------------------------------------------------------------
function s = tern_str(tf, a, b)
if tf, s = a; else, s = b; end
end

% -----------------------------------------------------------------------------
function s = pass(tf)
if tf, s = 'OK'; else, s = '<-- FAILS'; end
end
