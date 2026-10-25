function [y, info] = emg_filter(x, fs, varargin)
%EMG_FILTER  The canonical sEMG conditioning chain for MC 4206.
%
%   [y, info] = EMG_FILTER(x, fs)
%   [y, info] = EMG_FILTER(x, fs, 'band', [20 450], 'mains', 50, 'harmonics', 3)
%
%   x   [nsamples x nchannels] RAW EMG.  Run emg_quality FIRST.
%   y   filtered EMG, same size.
%   -----------------------------------------------------------------------------
%   THESE PARAMETERS ARE NOT ARBITRARY.  DO NOT CHANGE THEM WITHOUT A REASON
%   WRITTEN INTO THE DECISION LOG.
%   LITERATURE-VALIDATED 2026-09-25.  See DECISION_LOG D-25 for the full pass;
%   summary below.  
%   -----------------------------------------------------------------------------
%   BAND 20-450 Hz, 4th-order Butterworth, zero-phase (filtfilt).
%     Low cut 20 Hz  - De Luca, Gilmore, Kuznetsov & Roy (2010), "Filtering
%                      the surface EMG signal: Movement artifact and baseline
%                      noise contamination", J. Biomech. 43(8):1573-1579.
%                      They tested corner frequencies against movement
%                      artefact rejection on tibialis anterior and first
%                      dorsal interosseous and recommend a 20 Hz corner
%                      (Butterworth, 12 dB/oct = 2nd order per pass) as "the
%                      best compromise" - consistent with the older ISEK
%                      recommendation (Winter et al. 1980), and stricter than
%                      some competing conventions (5 Hz Merletti; 10 Hz
%                      J. Electromyogr. Kinesiol.).  This project runs 4th
%                      order per pass through filtfilt (~8th-order-equivalent
%                      magnitude roll-off), which is deliberately sharper than
%                      De Luca's minimum-sufficient 2nd order, not a departure
%                      from it.  Removes motion artefact, cable sway and the
%                      electrode half-cell drift, all of which live below
%                      ~20 Hz.  Going lower lets motion artefact into the
%                      envelope, where it looks exactly like a slow contraction.
%     High cut 450 Hz- within the general practical range widely reported for
%                      surface EMG conditioning (commonly quoted as roughly
%                      10-20 Hz to 450-500 Hz for the high/low corners
%                      respectively; e.g. Noraxon's EMG signal-processing
%                      guidance quotes 20-450 Hz directly). Essentially all
%                      surface EMG power is below this; above it you are
%                      amplifying amplifier noise.
%     Zero-phase     - filtfilt runs the filter forwards and backwards, so there
%                      is NO phase distortion and NO group delay.  That matters
%                      because we later compare EMG onset timing against the
%                      control log.  It is an OFFLINE operation: it uses future
%                      samples, so it cannot run in real time.  The real-time
%                      path uses a causal single-pass filter and accepts the
%                      delay - stated explicitly that difference in the thesis.
%
%   MAINS NOTCH at 50 Hz and harmonics.
%     Sri Lanka is 230 V / 50 Hz.  Note that 50 Hz sits INSIDE the 20-450 band
%     and inside the dominant 50-150 Hz sEMG energy, so it cannot be removed by
%     the bandpass.  Each notch is narrow (Q high) to take as little real signal
%     as possible.  If emg_quality reports mains_ratio above ~0.05 even after
%     notching, the fix is electrical so re-prep the skin, check the reference -
%     not more filtering.
%     Narrow notch + first few harmonics is standard engineering practice for
%     mains rejection in biosignals; this project did not find one single
%     canonical numeric standard for notch width, so the 2 Hz full width used
%     here is an engineering choice (tight enough to spare real signal at the
%     dominant 50-150 Hz sEMG energy) rather than a cited figure - stated
%     honestly as such per DECISION_LOG D-25.
%
%   ORDER OF OPERATIONS.  Band-pass, then notch, then (downstream) rectify and
%   envelope.  Never rectify before filtering: rectification is nonlinear and
%   folds high-frequency content down into the band you are about to keep.
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Revised on 2026-09-20 / Rev A


o.band      = [20 450];     % Hz   De Luca et al. 2010 / general practice (D-25)
o.order     = 4;            %      Butterworth order (per direction)
o.mains     = 50;           % Hz   [] to disable
o.harmonics = 3;            %      notch 50, 100, 150 Hz
o.notch_bw  = 2;            % Hz   full width of each notch
o.causal    = false;        %      true -> single-pass (real-time equivalent)
o = parseopts(o, varargin);

[N, nch] = size(x);        
nyq = fs/2;
hi  = min(o.band(2), 0.95*nyq);
lo  = max(o.band(1), 0.5);
if hi <= lo
    error('emg_filter:badBand','Band [%g %g] is invalid at fs = %g.', lo, hi, fs);
end

% ---- 1. band-pass ----------------------------------------------------------
[b, a] = butter(o.order, [lo hi]/nyq, 'bandpass');
y = applyf(b, a, x, o.causal);

% ---- 2. mains notches ------------------------------------------------------
notches = [];
if ~isempty(o.mains)
    for h = 1:o.harmonics
        f0 = o.mains*h;
        if f0 >= 0.95*nyq || f0 <= lo || f0 >= hi, continue; end
        wn = [f0-o.notch_bw/2, f0+o.notch_bw/2]/nyq;
        [bn, an] = butter(2, wn, 'stop');
        y = applyf(bn, an, y, o.causal);
        notches(end+1) = f0; 
    end
end

info = struct('band',[lo hi],'order',o.order,'notches',notches, ...
              'causal',o.causal,'fs',fs, ...
              'phase', tern(o.causal,'causal single-pass (has delay)', ...
                                     'zero-phase filtfilt (offline only)'));

if nargout == 0
    fprintf('\n  emg_filter\n');
    fprintf('  ------------------------------------------------------\n');
    fprintf('  band      %g - %g Hz, Butterworth order %d\n', lo, hi, o.order);
    fprintf('  notches   %s Hz\n', mat2str(notches));
    fprintf('  phase     %s\n', info.phase);
    fprintf('  ------------------------------------------------------\n\n');
end
end

% ==============================================================================
function y = applyf(b, a, x, causal)
if causal
    y = filter(b, a, x);
else
    y = filtfilt(b, a, x);
end
end

function s = tern(tf,a,b)
if tf, s = a; else, s = b; end
end

function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_filter:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
