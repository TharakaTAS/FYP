function B = emg_batch_mendeley(folder, varargin)
%EMG_BATCH_MENDELEY  Run the loader over several participants and CROSS-CHECK
%                    that the column-to-muscle mapping is consistent.
%
%   B = EMG_BATCH_MENDELEY('C:\Users\ashen\Downloads')
%   B = EMG_BATCH_MENDELEY(folder, 'files', {'10_raw.csv','11_raw.csv'})
%   B = EMG_BATCH_MENDELEY(folder, 'pattern', '*_raw.csv')
%
%   ---------------------------------------------------------------------------
%   WHY RUNNING ONE FILE IS NOT ENOUGH
%
%   A single participant tells you which column looked like which muscle IN THAT
%   RECORDING.  It cannot tell you whether the electrode order was the same for
%   every participant.  Those are different claims, and only the second one lets
%   you use the other 35 participants.
%
%     - Mapping the SAME in every file  -> the column order is a property of the
%       dataset's protocol.  You may rely on it, and you say so with evidence.
%     - Mapping FLIPS between files     -> the order is per-session.  Every file
%       must be identified individually, and any published result that assumed a
%       fixed order is wrong.
%
%   This is exactly the check that "just take columns 1 and 3" skips.
%
%   OUTPUT B
%     .file, .map, .conf, .aligned, .qpass, .corr_max, .mode_a   per file
%     .consistent    true if every file gave the same assignment
%     .table         summary, one row per file
%
%   MC 4206 / Team MC 04 / T. A. S. Tharaka / Rev A 2026-09-20

addpath(fileparts(mfilename('fullpath')));

o.files   = {};
o.pattern = '*_raw.csv';
o.plot    = false;
o.verbose = false;        % true -> show each file's full report
o = parseopts(o, varargin);

if isempty(o.files)
    d = dir(fullfile(folder, o.pattern));
    o.files = {d.name};
end
if isempty(o.files)
    error('emg_batch_mendeley:noFiles', ...
          'No files matching "%s" in %s', o.pattern, folder);
end
[~, ord] = sort(o.files);  o.files = o.files(ord);

nf = numel(o.files);
fprintf('\n=================================================================\n');
fprintf('  BATCH  -  %d files from %s\n', nf, folder);
fprintf('=================================================================\n');

B = struct('file',{},'map',{},'conf',{},'aligned',{},'qpass',{}, ...
           'corr_max',{},'mode_a',{},'ok',{});

for k = 1:nf
    f = fullfile(folder, o.files{k});
    fprintf('\n[%d/%d] %s ', k, nf, o.files{k});
    try
        if o.verbose
            D = emg_load_mendeley(f, 'plot', o.plot);
        else
            % Swallow the per-file report; the summary below is what matters.
            evalc('D = emg_load_mendeley(f, ''plot'', o.plot);');
        end
        fprintf('OK\n');
        B(end+1).file = o.files{k};                      %#ok<AGROW>
        B(end).map      = {D.map.muscle};
        B(end).conf     = D.conf;
        B(end).aligned  = D.aligned;
        B(end).qpass    = D.Q.allpass;
        c = D.corr; B(end).corr_max = max(c(~eye(size(c,1))));
        B(end).mode_a   = D.mode_a;
        B(end).ok       = isfield(D,'map_ok') && D.map_ok;
        B(end).anchored = D.n_anchored;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('        at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

if isempty(B)
    error('emg_batch_mendeley:allFailed','Every file failed to load.');
end

% ---- summary ---------------------------------------------------------------
fprintf('\n=================================================================\n');
fprintf('  COLUMN -> MUSCLE, PER FILE\n');
fprintf('=================================================================\n');
fprintf('  %-16s %-6s %-6s %-6s %-6s   %-7s %-7s %s\n', ...
        'file','col1','col2','col3','col4','aligned','maxcorr','verdict');
for k = 1:numel(B)
    fprintf('  %-16s %-6s %-6s %-6s %-6s   %-7s %7.2f %s\n', ...
            B(k).file, B(k).map{1}, B(k).map{2}, B(k).map{3}, B(k).map{4}, ...
            tern(B(k).aligned,'yes','NO'), B(k).corr_max, ...
            tern(B(k).ok,'PROCEED','check'));
end

fprintf('\n  CONFIDENCE PER COLUMN (%%)  -  >=35 is usable\n');
fprintf('  %-16s %7s %7s %7s %7s\n','file','col1','col2','col3','col4');
C = zeros(numel(B),4);
th_outlier = 20;     % mean confidence below this -> the WHOLE FILE is unreliable
outlier = false(numel(B),1);
for k = 1:numel(B)
    C(k,:) = 100*B(k).conf;
    outlier(k) = mean(C(k,:)) < th_outlier;
    fprintf('  %-16s', B(k).file);
    for j = 1:4
        fprintf('%6.0f%s', C(k,j), tern(C(k,j)>=35,' ','*'));
    end
    if outlier(k)
        fprintf('   <- OUTLIER (all 4 cols < %d%%)', th_outlier);
    end
    fprintf('\n');
end
fprintf('  %-16s', 'MEAN');
fprintf('%6.0f ', mean(C,1)); fprintf('\n');
fprintf('  (* = below 35 %%, do not treat that column as an independent muscle)\n');

% cycles-anchored breakdown, per file (this used to be a dead loop that
% printed nothing - kept as its own small table instead of trying to
% cram a single number into the column-per-muscle table above).
fprintf('\n  CYCLES ANCHORED (of 5 expected)\n');
for k = 1:numel(B)
    tag = '';
    if isfield(B(k),'anchored') && B(k).anchored < 5
        tag = '   <- fewer than 5, segmentation may be unreliable';
    end
    fprintf('  %-16s %d/5%s\n', B(k).file, B(k).anchored, tag);
end

% ---- consistency -----------------------------------------------------------
maps  = vertcat(B.map);
files = {B.file};
nf2   = numel(B);

% Per-column consistency.  An all-or-nothing test is useless here: if ONE
% column is ambiguous the whole verdict collapses, even when the columns you
% actually intend to use are identical in every file.  Judge each column,
% AND when a column is not unanimous, report the MAJORITY answer and name
% whoever dissents - a 4-of-5 split and a 2-2-1 split are different findings
% and must not print identically as "varies - not reliable".
colconst = false(1,4);
majority = cell(1,4);
majcount = zeros(1,4);
for j = 1:4
    u = unique(maps(:,j));
    colconst(j) = numel(u) == 1;
    counts = cellfun(@(v) sum(strcmp(maps(:,j), v)), u);
    [majcount(j), im] = max(counts);
    majority{j} = u{im};
end
consistent = all(colconst);
B(1).consistent = consistent; %#ok<STRNU>

% Column indices come from the MAJORITY map, not file 1's row - file 1 is
% not special, and using it meant an unlucky first file could point at the
% wrong column.
iE = find(strcmp(majority,'ECR'),1);
iF = find(strcmp(majority,'FCR'),1);

fprintf('\n=================================================================\n');
fprintf('  CONSISTENCY VERDICT\n');
fprintf('=================================================================\n');
if consistent
    fprintf('  The mapping is IDENTICAL across all %d files:\n', nf2);
    fprintf('    col1=%s  col2=%s  col3=%s  col4=%s\n', maps{1,:});
    fprintf('\n  >> The column order is a property of the dataset protocol, not\n');
    fprintf('     of the individual session.  You may apply it to the other\n');
    fprintf('     participants - and you can say so WITH EVIDENCE, which is a\n');
    fprintf('     different and much stronger claim than assuming it.\n');
    if ~isempty(iF) && ~isempty(iE)
        fprintf('\n     MODE A (D-13):  ECR = col %d,  FCR = col %d\n', iE, iF);
        fprintf('     two_ch = D.Xf(:, [%d %d]);\n', iE, iF);
    end
else
    fprintf('  Per-column result across %d files:\n', nf2);
    for j = 1:4
        if colconst(j)
            fprintf('     col%d : %-4s  CONSISTENT in every file\n', j, majority{j});
        elseif majcount(j) >= nf2 - 1
            dissenters = files(~strcmp(maps(:,j)', majority{j}));
            fprintf('     col%d : %-4s in %d/%d files; dissents: %s\n', ...
                    j, majority{j}, majcount(j), nf2, strjoin(dissenters, ', '));
        else
            u = unique(maps(:,j));
            fprintf('     col%d : varies (%s)  - not reliable\n', j, strjoin(u', ', '));
        end
    end

    modeA_strict   = ~isempty(iF) && ~isempty(iE) && colconst(iF) && colconst(iE);
    modeA_majority = ~isempty(iF) && ~isempty(iE) && colconst(iF) && majcount(iE) >= nf2-1;

    if modeA_strict
        fprintf('\n  >> The MODE A PAIR IS CONSISTENT even though other columns\n');
        fprintf('     are not.  ECR = col %d and FCR = col %d in every file.\n', iE, iF);
        fprintf('     You may use that pair, and you must state that the other\n');
        fprintf('     two columns could not be separated.\n');
        fprintf('\n     two_ch = D.Xf(:, [%d %d]);   %% [ECR FCR]\n', iE, iF);
    elseif modeA_majority
        dissenters = files(~strcmp(maps(:,iE)', majority{iE}));
        fprintf('\n  >> FCR = col %d is consistent in EVERY file.  ECR = col %d is\n', iF, iE);
        fprintf('     consistent in %d of %d files - the exception is %s, which\n', ...
                majcount(iE), nf2, strjoin(dissenters, ', '));
        fprintf('     must be identified individually (do not assume col %d there\n', iE);
        fprintf('     without checking).  For every other file:\n');
        fprintf('\n     two_ch = D.Xf(:, [%d %d]);   %% [ECR FCR]\n', iE, iF);
    else
        fprintf('\n  >> Electrode order is PER SESSION.  Identify every file\n');
        fprintf('     individually before use, and any result that assumed a\n');
        fprintf('     fixed column order is invalid.  This is a finding - record it.\n');
    end

    if any(outlier)
        fprintf('\n  >> Flagged OUTLIER above (all 4 columns below %d%% confidence): %s\n', ...
                th_outlier, strjoin(files(outlier), ', '));
        fprintf('     Its columns can land on the majority answer by chance at that\n');
        fprintf('     confidence level.  Treat it as unverified even where it agrees,\n');
        fprintf('     and confirm its channel identity by eye before using it.\n');
    end
end

nlow = sum(mean(C,1) < 35);
if nlow > 0
    fprintf('\n  NOTE: %d of 4 columns average below 35 %% confidence.  Expect\n', nlow);
    fprintf('  crosstalk between neighbouring electrodes - the max envelope\n');
    fprintf('  correlation across files is %.2f.  Report this as a limitation\n', ...
            max([B.corr_max]));
    fprintf('  of the dataset, and do not claim four independent muscles.\n');
end
fprintf('=================================================================\n\n');
end

% =============================================================================
function s = tern(tf,a,b)
if tf, s = a; else, s = b; end
end

function o = parseopts(o, args)
for k = 1:2:numel(args)
    if ~isfield(o, args{k})
        error('emg_batch_mendeley:badOption','Unknown option "%s".', args{k});
    end
    o.(args{k}) = args{k+1};
end
end
