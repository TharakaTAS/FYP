%RUN_ME_TODAY  Validate the whole processing pipeline against synthetic ground
%              truth - BEFORE any hardware exists.
%
%   Usage:  cd to D:\Tharaka\FYP\control, then just type
%
%       RUN_ME_TODAY
%
%   It needs no hardware, no subject, no data.  It takes about ten seconds and
%   produces THREE figures and a printed report.  (wrist_params prints a table
%   but draws nothing - that is why it is three and not four.)
%
%   ---------------------------------------------------------------------------
%   WHY THIS IS NOT A TOY
%
%   The surrogate-first principle, applied to software.  When a real
%   identification returns k_theta = 95 N.m/rad, you need to already know that
%   the code recovers 95 when the true answer IS 95.  Otherwise you cannot tell
%   a wrist from a bug.
%
%   So: generate data with a KNOWN answer, run the real analysis scripts on it
%   unmodified, and check that what comes out is what went in.  That is a result
%   you can show on Wednesday, and it is honest - it validates the method, not
%   the device.
%
%   WHAT IT PROVES
%     1. The parameter set passes its own design-rule checks
%     2. wrist_id_ktheta.m recovers a known transmission stiffness
%     3. wrist_id_endfeel.m recovers a known end-feel index
%     4. emg_baseline.m turns two EMG channels into a signed intent torque
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev D 2026-09-20

clear; close all; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'identification'));
rng(4206);                              % reproducible

fprintf('\n');
fprintf('===========================================================\n');
fprintf('  MC 4206  -  PIPELINE VALIDATION AGAINST SYNTHETIC TRUTH\n');
fprintf('  No hardware required.  Run date: %s\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('===========================================================\n');

%% ------------------------------------------------------------------ 1. PARAMS
fprintf('\n[1/4]  Parameter set and design-rule checks\n');
wrist_params('normal');                 % prints its own check table

P = wrist_params('normal');

%% ------------------------------------------------- 2. TEST 0  (k_theta)
fprintf('\n[2/4]  Test 0  -  can we recover a KNOWN cable stiffness?\n\n');

k_true    = 95;          % N.m/rad   <- the answer the script must find
hyst      = 0.06;        % 6% difference between up and down sweeps
lever     = 0.100;       % m
r_m       = 0.020;  r_j = 0.040;
q         = 360/4096;    % encoder LSB, degrees

% MASSES.  The first version of this demo used 50-200 g and FAILED, which is how
% we found the defect now recorded as F-08: at k = 95, 200 g on a 100 mm lever
% gives 0.118 deg of full-scale wind-up = 1.3 encoder counts.  You cannot fit a
% line through 1.3 counts.  1 kg gives 0.59 deg = 6.7 counts.
masses    = [0 200 400 600 800 1000];

% DITHER.  A perfectly still rig quantises deterministically and averaging does
% NOT help.  Real rigs have mechanical noise, and the protocol now adds a small
% deliberate servo oscillation during the averaging window.  Both are modelled
% here as noise of about 1.5 LSB, which is what makes the 5 s average resolve
% below one count.
dither_deg = 1.5*q;
navg       = 2500;       % 5 s at 500 Hz

rows = {};
for dr = [1 -1]
    k_dir = k_true * (1 + dr*hyst/2);
    for m = masses
        tau   = (m/1000)*9.80665*lever*dr;        % N.m, from first principles
        Delta = tau / k_dir;                        % rad of wind-up
        th_m  = 30;                                 % deg, servo holds position
        th_mr = th_m * (r_m/r_j);                   % deg, referred to joint
        th_j  = th_mr - rad2deg(Delta);
        % Quantise WITH dither, then average - this is what the real protocol
        % does over a 5 s window at 500 Hz.
        th_m_q = mean(round((th_m + dither_deg*randn(navg,1))/q)*q);
        th_j_q = mean(round((th_j + dither_deg*randn(navg,1))/q)*q);
        rows(end+1,:) = {m, dr, lever, th_m_q, th_j_q, 0.01, 0.01, ''}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', ...
      {'mass_g','direction','lever_m','theta_m_deg','theta_j_deg', ...
       'theta_m_std','theta_j_std','notes'});

csv = fullfile(tempdir, 'SYNTH_test0_ktheta.csv');
writetable(T, csv);

R0 = wrist_id_ktheta(csv, 'r_m', r_m, 'r_j', r_j, 'plot', true);

err0 = 100*(R0.k_theta - k_true)/k_true;
fprintf('  >>> TRUTH  k_theta = %.1f N.m/rad\n', k_true);
fprintf('  >>> FOUND  k_theta = %.1f N.m/rad   ERROR %+.2f %%   %s\n\n', ...
        R0.k_theta, err0, verdict(abs(err0) < 3));

%% ------------------------------------------------- 3. TEST 2b (end-feel)
fprintf('[3/4]  Test 2b  -  can we recover a KNOWN end-feel?\n\n');

K_mid_true  = 0.90;                 % N.m/rad, linear region
th_knee     = deg2rad(18);          % rad, where it starts to harden
stiffen     = 22;                   % N.m/rad^2, quadratic hardening

th   = deg2rad(linspace(0, 32, 600))';
tauQ = K_mid_true*th + stiffen*max(th - th_knee, 0).^2;
tauQ = tauQ + 0.004*randn(size(tauQ));            % torque-estimate noise

R2b = wrist_id_endfeel(rad2deg(th), tauQ, 'plot', true);

fprintf('  >>> TRUTH  K_mid = %.2f N.m/rad, knee at %.0f deg\n', ...
        K_mid_true, rad2deg(th_knee));
fprintf('  >>> FOUND  K_mid = %.2f N.m/rad, knee at %.0f deg, index %.2f  %s\n\n', ...
        R2b.K_mid, R2b.theta_knee_deg, R2b.endfeel_idx, ...
        verdict(abs(R2b.K_mid-K_mid_true)/K_mid_true < 0.15));

%% ------------------------------------------------- 4. EMG BASELINE
fprintf('[4/4]  EMG baseline  -  two channels in, signed intent torque out\n\n');

fsE = 1000;  tE = (0:1/fsE:12-1/fsE)';
burst = @(t0,t1,amp) amp*((tE>=t0)&(tE<t1));

% agonist bursts at 3-5 s and 8-10 s; antagonist burst at 6-7 s (co-contraction)
drive_ago = burst(3,5,1.0) + burst(8,10,0.55) + burst(6,7,0.8);
drive_ant = burst(6,7,0.8);

env2raw = @(d) (0.02 + 0.55*d) .* randn(size(tE));   % amplitude-modulated noise
emg_ago = env2raw(drive_ago) + 0.004*sin(2*pi*50*tE);   % + mains
emg_ant = env2raw(drive_ant) + 0.004*sin(2*pi*50*tE);

% TWO runs, and the comparison is the point.
%
%   (a) auto-normalised: each channel divided by its OWN 95th percentile.  The
%       agonist channel has larger bursts elsewhere in the record, so its
%       reference is larger, so during an EQUAL co-contraction the antagonist
%       looks relatively stronger and the intent goes NEGATIVE.  That is not a
%       bug in the algorithm - it is what happens without a proper MVC.
%
%   (b) MVC-normalised: both channels divided by a measured maximum.  Equal
%       drive now cancels, and the intent correctly sits at zero.
%
% The lesson, which goes in the protocol: RECORD A PER-CHANNEL MVC (or a
% standard calibration contraction) AT THE START OF EVERY SESSION.
[tauA, infoA] = emg_baseline(emg_ago, emg_ant, fsE, 'rest_win', [1 2*fsE]);

mvc = 0.02 + 0.55*1.0;                       % the true full-drive amplitude
[tauE, infoE] = emg_baseline(emg_ago, emg_ant, fsE, ...
                             'mvc_ago', 0.62*mvc, 'mvc_ant', 0.62*mvc, ...
                             'rest_win', [1 2*fsE]);

figure('Name','EMG baseline - proportional myoelectric intent','Color','w');
subplot(3,1,1); plot(tE, emg_ago, 'LineWidth',0.4); grid on; box on;
  ylabel('agonist [V]'); title('Raw EMG (synthetic)');
subplot(3,1,2); hold on; grid on; box on;
  plot(tE, infoE.a_ago, 'LineWidth',1.2, 'DisplayName','agonist activation');
  plot(tE, infoE.a_ant, 'LineWidth',1.2, 'DisplayName','antagonist activation');
  ylabel('normalised'); legend('Location','northwest');
subplot(3,1,3); hold on; grid on; box on;
  plot(tE, tauA, 'LineWidth',1.2, 'DisplayName','auto-normalised (no MVC)');
  plot(tE, tauE, 'LineWidth',1.6, 'DisplayName','MVC-normalised');
  yline(0,'-','HandleVisibility','off');
  ylabel('\tau_{intent} [N.m]'); xlabel('time [s]'); legend('Location','northwest');
  title('Intent torque - co-contraction at 6-7 s cancels ONLY with a proper MVC');

coA = mean(abs(tauA(tE>=6.2 & tE<=6.8)));
co  = mean(abs(tauE(tE>=6.2 & tE<=6.8)));
ag  = mean(tauE(tE>=3.5 & tE<=4.5));
lim = 0.15*max(abs(tauE));
fprintf('  >>> agonist-only burst         ->  mean intent %+.3f N.m\n', ag);
fprintf('  >>> co-contraction, auto-norm  ->  mean |intent| %.3f N.m   %s\n', coA, ...
        verdict(coA < lim));
fprintf('  >>> co-contraction, MVC-norm   ->  mean |intent| %.3f N.m   %s\n', co, ...
        verdict(co < lim));
fprintf('      A braced limb is NOT a request to move.  The antagonist channel\n');
fprintf('      cancels it - but only if BOTH channels are normalised against a\n');
fprintf('      real maximum.  Auto-normalising by a percentile of the record\n');
fprintf('      gives each channel a different reference and the cancellation\n');
fprintf('      fails.  PROTOCOL: record a per-channel MVC every session.\n\n');

%% ------------------------------------------------------------------ SUMMARY
fprintf('===========================================================\n');
fprintf('  SUMMARY  -  what this run demonstrates\n');
fprintf('===========================================================\n');
fprintf('  Test 0  k_theta recovery ................ %+.2f %% error  %s\n', ...
        err0, verdict(abs(err0) < 3));
fprintf('  Test 2b K_mid recovery ................. %+.2f %% error  %s\n', ...
        100*(R2b.K_mid-K_mid_true)/K_mid_true, ...
        verdict(abs(R2b.K_mid-K_mid_true)/K_mid_true < 0.15));
fprintf('  Test 2b end-feel index ................. %.2f  [%s]\n', ...
        R2b.endfeel_idx, R2b.class);
fprintf('  EMG co-contraction, auto-normalised .... %s  (expected)\n', ...
        verdict(coA < lim));
fprintf('  EMG co-contraction, MVC-normalised ..... %s\n', ...
        verdict(co < lim));
fprintf('-----------------------------------------------------------\n');
fprintf('  The analysis pipeline recovers known truth before any\n');
fprintf('  hardware exists.  When real data arrives, a wrong answer\n');
fprintf('  means the WRIST is interesting - not that the code is broken.\n');
fprintf('\n  Two defects were found by THIS demo, before the bench:\n');
fprintf('    F-08  Test 0 masses were 5x too light - 1.3 encoder counts of\n');
fprintf('          full-scale wind-up.  Now 1 kg, plus deliberate dither.\n');
fprintf('    F-09  EMG auto-normalisation breaks co-contraction rejection.\n');
fprintf('          A per-channel MVC is now a protocol requirement.\n');
fprintf('===========================================================\n\n');


function s = verdict(tf)
if tf, s = 'PASS'; else, s = '<-- FAIL'; end
end
