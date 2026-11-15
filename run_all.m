%% run_all.m - make the figures and tables of the article
% Salivary and plasma cytokine kinetics during experimental human endotoxaemia:
% macrophage inflammatory protein-1alpha and -1beta as candidate non-invasive biomarkers.
%
% Put the two data files in data/ (see data/README.md) and run this file.
% Output: results/figures (Figures 2-6 and S2-S8) and results/tables (see README.md).

REFIT  = false;     % false: use the fit of the article (fit/endotox_fit.mat), about 15-30 min
                    % true:  refit the model with 200 bootstrap refits (about 3 h with 16 workers)
PYTHON = "python3"; if ispc, PYTHON = "python"; end      % or the full path to Python

cd(fileparts(mfilename('fullpath')));
DATA   = fullfile('data', 'DEDICATE-LPS SALIVA.xlsx');
FIGDIR = fullfile('results', 'figures');
TABDIR = fullfile('results', 'tables');
assert(isfile(DATA), 'Data file not found: put the data files in data/ (see data/README.md).');
if ~isfolder(TABDIR), mkdir(TABDIR); end

% 1. model fit (Tables S2-S4)
if REFIT
    [F, Tp, Ts, Tprot, C] = endotoxemia_model(DATA, 'save', fullfile('results', 'endotox_fit.mat'));
else
    [F, Tp, Ts, Tprot, C] = endotoxemia_model(DATA, 'fit', fullfile('fit', 'endotox_fit.mat'));
end
writetable(F.tables.model_metrics, fullfile(TABDIR, 'model_metrics.csv'));
writetable(F.tables.parameters, fullfile(TABDIR, 'parameters.csv'));
writetable(F.tables.baselines, fullfile(TABDIR, 'baselines.csv'));

% 2. figures and data tables
figures_tables

% 3. lag with each treatment group left out of the fit (Table S7; no bootstrap)
lags = @(Fx, label) [table(repmat(label, height(Fx.tables.lags), 1), 'VariableNames', {'fit'}), ...
    renamevars(Fx.tables.lags, {'Cytokine', 'lag', 'lagCentroid'}, {'cytokine', 'peakLag', 'centroidLag'})];
LOGO = lags(F, "all three arms");
pooled = ["placebo", "anakinra", "tocilizumab"];
for g = pooled
    Fg = endotoxemia_model(DATA, 'groups', pooled(pooled ~= g), 'nBoot', 0);
    LOGO = [LOGO; lags(Fg, "without " + g)]; %#ok<AGROW>
end
writetable(LOGO, fullfile(TABDIR, 'lags_leave_one_group_out.csv'));

% 4. Table 1, values below the LLOQ, Tables S12-S13
assert(system(PYTHON + " " + fullfile('python', 'tables.py')) == 0, ...
    'The Python step failed: install python/requirements.txt or set PYTHON above.');
fprintf('\nDone. Results are in %s\n', fullfile(pwd, 'results'));
