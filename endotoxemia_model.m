function [F, Tp, Ts, Tprot, C] = endotoxemia_model(dataFile, varargin)
%ENDOTOXEMIA_MODEL  Data import and compartmental plasma-to-saliva model of the cytokine
% response to LPS (10 cytokines: TNF, IL-6, IL-8, IL-10, IL-1RA, G-CSF, IP-10, MCP-1,
% MIP-1alpha, MIP-1beta).
%
%   [F, Tp, Ts, Tprot, C] = endotoxemia_model(dataFile)              fit + 200 bootstrap refits
%   ... = endotoxemia_model(dataFile, 'fit', 'fit/endotox_fit.mat')   load a saved fit
%   ... = endotoxemia_model(dataFile, 'groups', ["placebo","anakinra"], 'nBoot', 0)
%   ... = endotoxemia_model(dataFile, 'save', file)                   save the new fit
%
% Output
%   F       the fit: parameters (F.p, F.qJoint), simulated curves (F.sim), derived metrics
%           (F.metrics, columns as metric_names), bootstrap results (F.unc) and the model
%           tables (F.tables: metrics, lags, parameters, baselines)
%   Tp, Ts  cleaned plasma and saliva data (id, time in min, cyto, conc in pg/mL, group)
%   Tprot   salivary total protein
%   C       model settings
%
% The fit has three stages: (A) the shared upstream kinetics with TNF, (B) every other
% cytokine with the upstream part fixed, (C) all parameters jointly. Uncertainty: the
% participants are resampled with replacement within each group and stage C is refitted.
% Requires the Optimization Toolbox; the Parallel Computing Toolbox is optional.

opt = struct('groups', ["placebo", "anakinra", "tocilizumab"], 'nBoot', 200, ...
             'fit', '', 'save', '', 'parallel', true);
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i + 1}; end
SEED = 2026; CI_LEVEL = 0.95; BOOT_MAXITER = 150;

[Tp, Ts, Tprot] = load_data(dataFile);

C = struct();
C.t_end = 1440;                                     % minutes (24 h)
C.tGrid = linspace(0, C.t_end, 5761);
C.cytokineNames = ["TNF", "IL6", "IL8", "IL10", "IL1RA", "GCSF", "IP10", "MCP1", "MIP1a", "MIP1B"];
C.hasLocalFun = @(nm) strcmp(char(nm), 'IL8');      % local oral production term: IL-8 only
% the drug acts on the cytokine itself: IL-1RA without anakinra, IL-6 without tocilizumab
C.exclusionGroupMap = containers.Map({'IL1RA', 'IL6'}, {'anakinra', 'tocilizumab'});
% time windows (min) that are up-weighted in the fit
C.windowMap = struct('TNF', [30 90], 'IL6', [120 210], 'IL8', [120 210], 'IL10', [120 210], ...
    'IL1RA', [60 150], 'GCSF', [120 210], 'IP10', [180 300], 'MCP1', [90 180], ...
    'MIP1a', [60 150], 'MIP1B', [60 150]);
C.peakWeight = 4;                                   % weight of the peak-height residuals
C.midWindowWeight = 3;                              % weight of the points in windowMap
C.stageCMaxEval = 30000;
C.useParallel = false;
C.lloqS = struct();                                 % salivary LLOQ per cytokine (figures)
for nm = C.cytokineNames
    v = Ts.lloq(Ts.cyto == nm); C.lloqS.(char(nm)) = median(v(isfinite(v)));
end

unit = struct('label', strjoin(opt.groups, " + "), 'groups', opt.groups);
if ~isempty(opt.fit)
    S = load(opt.fit);
    if isfield(S, 'F'), F = S.F; else, F = S.fits{1}; end
    fprintf('Loaded the fit from %s\n', opt.fit);
else
    C.useParallel = opt.parallel && start_pool();
    rng(SEED);
    fprintf('\nFitting: %s\n', unit.label);
    F = fit_unit(Tp, Ts, unit, C);
    if opt.nBoot > 0
        F.unc = uncertainty_bootstrap(Tp, Ts, unit, F, C, opt.nBoot, BOOT_MAXITER, CI_LEVEL, SEED + 1);
    end
    if ~isempty(opt.save), save(opt.save, 'F', '-v7.3'); end
end
F.tables = model_tables(F, CI_LEVEL);
end


%% ========================================================================
%  DATA
%  ========================================================================
function [Tp, Ts, Tprot] = load_data(file)
% Only the study day is analysed: Day 0 and its 24-h sample (coded D0 and D1).
DAYS = ["D0", "D1"];

% ---- saliva (Olink export)
Ts = readtable(file, detectImportOptions(file, 'Sheet', "Olink flex saliva", 'VariableNamingRule', 'preserve'));
Ts = Ts(string(Ts.SampleType) == "SAMPLE", :);
qv = num(Ts.QuantifiedValue); ll = num(Ts.LLOQ);
% failed wells: >= 8 of the 10 cytokines below the LLOQ or not quantified
[sid, ~, g] = unique(string(Ts.SampleID));
failed = sid(accumarray(g, double(isnan(qv) | qv < ll)) >= 8);
keep = ~ismember(string(Ts.SampleID), failed);
Ts = Ts(keep, :); qv = qv(keep); ll = ll(keep);
fprintf('Saliva: excluded %d failed well(s): %s\n', numel(failed), strjoin(failed, ', '));
[id, day, tm] = parse_sample_id(Ts.SampleID);
assay = containers.Map({'TNF', 'IL6', 'CXCL8', 'IL10', 'IL1RN', 'CSF3', 'CXCL10', 'CCL2', 'CCL3', 'CCL4'}, ...
    {'TNF', 'IL6', 'IL8', 'IL10', 'IL1RA', 'GCSF', 'IP10', 'MCP1', 'MIP1a', 'MIP1B'});
known = isKey(assay, cellstr(string(Ts.Assay)));
cyto = strings(height(Ts), 1);
cyto(known) = string(values(assay, cellstr(string(Ts.Assay(known)))));
Ts = table(id, tm, cyto, qv, ll, 'VariableNames', {'id', 'time', 'cyto', 'conc', 'lloq'});
Ts = Ts(ismember(day, DAYS) & known & ~isnan(tm) & ~isnan(Ts.conc), :);
Ts.blq = ~isnan(Ts.lloq) & Ts.conc < Ts.lloq;       % below the LLOQ: LLOQ/sqrt(2)
Ts.conc(Ts.blq) = Ts.lloq(Ts.blq) / sqrt(2);

% ---- plasma (columns id, day, time, cyto, pg/mL, Group)
Tp = readtable(file, detectImportOptions(file, 'Sheet', "plasma cytokines ", 'VariableNamingRule', 'preserve'));
Tp.Properties.VariableNames = {'id', 'day', 'time', 'cyto', 'conc', 'group'};
Tp = table(string(Tp.id), string(Tp.day), num(Tp.time), string(Tp.cyto), num(Tp.conc), ...
    lower(strtrim(string(Tp.group))), 'VariableNames', {'id', 'day', 'time', 'cyto', 'conc', 'group'});
Tp.time(Tp.day == "D1") = 1440;
Tp = Tp(ismember(Tp.day, DAYS) & ~isnan(Tp.conc) & ~isnan(Tp.time), {'id', 'time', 'cyto', 'conc', 'group'});

% treatment group of the saliva participants
[u, ia] = unique(Tp.id, 'stable');
[tf, loc] = ismember(Ts.id, u);
Ts.group = strings(height(Ts), 1); Ts.group(tf) = Tp.group(ia(loc(tf)));

% outliers; salivary values below the LLOQ are never removed
keep = outlier_keep(Tp, false(height(Tp), 1));
fprintf('Plasma: removed %d of %d values as outliers\n', sum(~keep), height(Tp)); Tp = Tp(keep, :);
keep = outlier_keep(Ts, Ts.blq);
fprintf('Saliva: removed %d of %d values as outliers\n', sum(~keep), height(Ts)); Ts = Ts(keep, :);

% ---- salivary total protein (|z| > 3 per time point removed)
Tprot = readtable(file, detectImportOptions(file, 'Sheet', "saliva protein", 'VariableNamingRule', 'preserve'));
Tprot.Properties.VariableNames = {'sampleID', 'protein'};
Tprot.sampleID = string(Tprot.sampleID);
Tprot = Tprot(Tprot.sampleID ~= "" & ~isnan(Tprot.protein), :);
[id, day, tm] = parse_sample_id(Tprot.sampleID);
Tprot.id = id; Tprot.time = tm;
Tprot = Tprot(ismember(day, DAYS) & ~isnan(tm), :);
G = findgroups(Tprot.time); keep = true(height(Tprot), 1);
for gi = 1:max(G)
    idx = find(G == gi); v = Tprot.protein(idx);
    if numel(v) >= 3 && std(v) > 0, keep(idx(abs((v - mean(v)) / std(v)) > 3)) = false; end
end
Tprot = Tprot(keep, :);
end

function [id, day, tm] = parse_sample_id(s)
% "R34.72 D0_T3" -> id R34.72, day D0, 180 min; "R34.72 D1" -> 24 h (1440 min); else NaN
s = strrep(string(s), ' ', '_');
id = extractBefore(s, '_'); rest = extractAfter(s, '_');
day = extractBefore(rest + "_", '_');
tm = 60 * str2double(extractAfter(rest, "_T"));
tm(~contains(rest, '_')) = 1440;
end

function x = num(x)
if ~isnumeric(x), x = str2double(strrep(string(x), ',', '.')); end
end

function keep = outlier_keep(T, protect)
% Within cytokine x time x group: drop values whose log10 lies more than 3.5 robust SDs
% (1.4826 * MAD) from the median. Rows in protect are kept.
keep = true(height(T), 1);
G = findgroups(T.cyto, T.time, T.group);
for gi = 1:max(G)
    idx = find(G == gi); vals = T.conc(idx); ok = vals > 0;
    if sum(ok) < 4, continue; end
    lv = log10(vals(ok)); io = idx(ok);
    sd = 1.4826 * median(abs(lv - median(lv)));
    if sd > 0, keep(io(abs(lv - median(lv)) > 3.5 * sd)) = false; end
end
keep = keep | protect;
end


%% ========================================================================
%  MODEL EQUATIONS AND RESIDUALS
%  ========================================================================
function Hv = Hup(X, K, n)
Hv = X.^n ./ (K.^n + X.^n + eps);
end

function Hv = Hdown(X, K, n)
Hv = K.^n ./ (K.^n + X.^n + eps);
end

function dY = coreODE(~, Y, kLPS, kM, kM1)
dY = [-kLPS*Y(1); Y(1) - kM*Y(2); kM*Y(2) - kM1*Y(3)];
end

function dY = stageA_ODE(~, Y, q)
% q = [kLPS kM kM1 KLPS nLPS s kmRNA kC kBS kQ kS] (TNF)
dY = zeros(7, 1);
dY(1:3) = coreODE([], Y(1:3), q(1), q(2), q(3));
dY(4) = Hup(Y(1), q(4), q(5)) - q(7)*Y(4);          % mRNA
dY(5) = q(6)*Y(4)*Y(3) - q(8)*Y(5);                 % plasma
dY(6) = q(9)*Y(5) - q(10)*Y(6);                     % oral transit
dY(7) = q(10)*Y(6) - q(11)*Y(7);                    % saliva
end

function dY = cytokineODE(t, Y, q, KLPS, nLPS, LPSinterp, M1interp)
% q = [s kmRNA kC tDelay nDelay kBS kQ kS local localDelay]
M1now = M1interp(t);
dmrna = Hup(LPSinterp(t), KLPS, nLPS) - q(2)*Y(1);
dX = q(1)*Y(1)*M1now*Hup(t, q(4), q(5)) - q(3)*Y(2);
localX = q(9) * Hup(t, q(10), q(5)) * Hup(M1now, 0.2, 1);
dY = [dmrna; dX; q(6)*Y(2) + localX - q(7)*Y(3); q(7)*Y(3) - q(8)*Y(4)];
end

function dxdt = modelODE(t, Y, p, cytokineNames)
% States: LPS, M, M1, then per cytokine mRNA, plasma, oral transit, saliva (above baseline).
% TNF has tDelay = 0 and nDelay = 1, so its activation gate is 1 for t > 0.
nC = numel(cytokineNames);
dxdt = zeros(3 + 4*nC, 1);
dxdt(1:3) = coreODE([], Y(1:3), p.kLPS, p.kM, p.kM1);
M1 = Y(3);
LPSsignal = Hup(Y(1), p.KLPS, p.nLPS);
for k = 1:nC
    nm = char(cytokineNames(k)); cfg = p.(nm); i0 = 3 + (k-1)*4;
    X = Y(i0+2);
    if strcmp(nm, 'IL10')
        dmrna = LPSsignal * Hdown(X, 150, 2) - cfg.kmRNA * Y(i0+1);   % weak negative feedback
    else
        dmrna = LPSsignal - cfg.kmRNA * Y(i0+1);
    end
    localX = cfg.local * Hup(t, cfg.localDelay, cfg.nDelay) * Hup(M1, 0.2, 1);
    dxdt(i0+1) = dmrna;
    dxdt(i0+2) = cfg.s * Y(i0+1) * M1 * Hup(t, cfg.tDelay, cfg.nDelay) - cfg.kC * X;
    dxdt(i0+3) = cfg.kBS * X + localX - cfg.kQ * Y(i0+3);
    dxdt(i0+4) = cfg.kQ * Y(i0+3) - cfg.kS * Y(i0+4);
end
end

function r = fit_resid(model, data, w)
% Weighted log residuals (targets are geometric means).
r = sqrt(w(:)) .* (log(max(model(:), 1e-9)) - log(max(data(:), 1e-9)));
end

function ts = solve_grid(tspan, varargin)
% 2-min grid plus every sampling time.
tAll = vertcat(varargin{:});
tAll = tAll(tAll >= tspan(1) & tAll <= tspan(end));
ts = unique([(tspan(1):2:tspan(end))'; tspan(end); tAll(:)]);
end

function w = pointWeights(tVec, window, weight)
w = ones(size(tVec));
w(tVec >= window(1) & tVec <= window(2)) = weight;
end

function r = pair_resid(tSol, yP, yS, d, bP, bS, window, C)
% Residuals of one cytokine: plasma and saliva points, and both peak heights.
plasmaModel = bP + interp1(tSol, yP, d.tPlasma, 'linear', 'extrap');
salivaModel = bS + interp1(tSol, yS, d.tSaliva, 'linear', 'extrap');
r = [fit_resid(plasmaModel, d.mPlasma, pointWeights(d.tPlasma, window, C.midWindowWeight));
     fit_resid(salivaModel, d.mSaliva, pointWeights(d.tSaliva, window, C.midWindowWeight));
     fit_resid(bP + max(yP), max(d.mPlasma), C.peakWeight);
     fit_resid(bS + max(yS), max(d.mSaliva), C.peakWeight)];
end

function r = stageA_residuals(q, x0, odeOpts, d, bP, C)
[tS, Y] = ode45(@(t, Y) stageA_ODE(t, Y, q), solve_grid([0 C.t_end], d.tPlasma, d.tSaliva), x0, odeOpts);
r = pair_resid(tS, Y(:, 5), Y(:, 7), d, bP, d.baseline, C.windowMap.TNF, C);
end

function r = cytokine_residuals(q, odeOpts, KLPS, nLPS, LPSinterp, M1interp, d, bP, window, C)
[tS, Y] = ode45(@(t, Y) cytokineODE(t, Y, q, KLPS, nLPS, LPSinterp, M1interp), ...
    solve_grid([0 C.t_end], d.tPlasma, d.tSaliva), [0; 0; 0; 0], odeOpts);
r = pair_resid(tS, Y(:, 2), Y(:, 4), d, bP, d.baseline, window, C);
end

function [startIdx, totalLen] = cytokineParamLayout(cytokineNames, hasLocalFun)
% Stage-C vector: 5 upstream parameters, then 8 (10 with a local term) per cytokine.
startIdx = struct(); ptr = 6;
for k = 1:numel(cytokineNames)
    startIdx.(char(cytokineNames(k))) = ptr;
    ptr = ptr + 8 + 2*hasLocalFun(cytokineNames(k));
end
totalLen = ptr - 1;
end

function p = unpackJointParams(q, p, cytokineNames, hasLocalFun)
p.kLPS = q(1); p.kM = q(2); p.kM1 = q(3); p.KLPS = q(4); p.nLPS = q(5);
startIdx = cytokineParamLayout(cytokineNames, hasLocalFun);
for k = 1:numel(cytokineNames)
    nm = char(cytokineNames(k)); i = startIdx.(nm);
    p.(nm).s = q(i); p.(nm).kmRNA = q(i+1); p.(nm).kC = q(i+2);
    p.(nm).tDelay = q(i+3); p.(nm).nDelay = q(i+4);
    p.(nm).kBS = q(i+5); p.(nm).kQ = q(i+6); p.(nm).kS = q(i+7);
    if hasLocalFun(cytokineNames(k)), p.(nm).local = q(i+8); p.(nm).localDelay = q(i+9); end
end
end

function r = joint_residuals(q, p, x0, odeOpts, target, cytokineNames, C)
p = unpackJointParams(q, p, cytokineNames, C.hasLocalFun);
tData = [];
for k = 1:numel(cytokineNames)
    d = target.(char(cytokineNames(k))); tData = [tData; d.tPlasma(:); d.tSaliva(:)]; %#ok<AGROW>
end
[tS, Y] = ode45(@(t, Y) modelODE(t, Y, p, cytokineNames), solve_grid([0 C.t_end], tData), x0, odeOpts);
r = [];
for k = 1:numel(cytokineNames)
    nm = char(cytokineNames(k)); i0 = 3 + (k-1)*4;
    r = [r; pair_resid(tS, Y(:, i0+2), Y(:, i0+4), target.(nm), p.plasmaBaseline.(nm), ...
        p.salivaBaseline.(nm), C.windowMap.(nm), C)]; %#ok<AGROW>
end
end


%% ========================================================================
%  FITTING
%  ========================================================================
function R = get_rows(T, groups, cyto, C)
% Rows of one cytokine in the given groups, with the cytokine-specific exclusion.
R = T(ismember(T.group, groups) & T.cyto == cyto, :);
if isKey(C.exclusionGroupMap, char(cyto))
    R = R(R.group ~= string(C.exclusionGroupMap(char(cyto))), :);
end
end

function [target, cytoKeep] = build_targets(Tp, Ts, unit, cytoNames, C)
% Fit targets: geometric mean per sampling time, 0-24 h. Baseline: the first time point
% after LPS (0.5 h) of each matrix.
gm = @(x) exp(mean(log(x), 1));
target = struct(); cytoKeep = strings(1, 0);
for nm = cytoNames
    Cp = get_rows(Tp, unit.groups, nm, C); Cp = Cp(Cp.conc > 0 & Cp.time >= 0 & Cp.time <= C.t_end, :);
    Cs = get_rows(Ts, unit.groups, nm, C); Cs = Cs(Cs.conc > 0 & Cs.time >= 0 & Cs.time <= C.t_end, :);
    if isempty(Cp) || isempty(Cs), continue; end
    [G, tp] = findgroups(Cp.time); mp = splitapply(gm, Cp.conc, G);
    [G, ts] = findgroups(Cs.time); ms = splitapply(gm, Cs.conc, G);
    target.(char(nm)) = struct('tPlasma', tp, 'mPlasma', mp, 'tSaliva', ts, 'mSaliva', ms, ...
        'baseline', ms(1), 'baselinePlasma', mp(1));
    cytoKeep(end+1) = nm; %#ok<AGROW>
end
end

function p = initial_params(target, cytoNames, C)
p = struct('LPS0', 1.0, 'kLPS', 0.030, 'kM', 0.050, 'kM1', 0.020, 'KLPS', 0.05, 'nLPS', 1.0);
for nm = C.cytokineNames
    p.(char(nm)) = struct('s', 0.5, 'kmRNA', 0.20, 'kC', 0.02, 'tDelay', 60, 'nDelay', 4, ...
        'kBS', 0.001, 'kQ', 0.05, 'kS', 0.15, 'local', 0, 'localDelay', 0);
end
p.TNF = struct('s', 1.00, 'kmRNA', 0.30, 'kC', 0.030, 'tDelay', 0, 'nDelay', 1, ...
    'kBS', 0.00065, 'kQ', 0.050, 'kS', 0.12, 'local', 0, 'localDelay', 0);
p.IL6.s = 0.69; p.IL6.kC = 0.008; p.IL6.tDelay = 50; p.IL6.kBS = 0.00150; p.IL6.kQ = 0.08; p.IL6.kS = 0.18;
p.IL8.s = 0.133; p.IL8.kC = 0.008; p.IL8.kBS = 0.040; p.IL8.kQ = 0.020; p.IL8.kS = 0.12;
p.IL8.local = 180; p.IL8.localDelay = 90;
p.IL10.s = 0.56; p.IL10.kmRNA = 0.18; p.IL10.kC = 0.008; p.IL10.tDelay = 50;
p.IL10.kBS = 0.00040; p.IL10.kQ = 0.10; p.IL10.kS = 0.18;
p.IL1RA.kC = 0.010; p.IL1RA.tDelay = 40;  p.GCSF.kC = 0.006;  p.GCSF.tDelay = 80;
p.IP10.kC  = 0.004; p.IP10.tDelay  = 100; p.MCP1.kC = 0.010;  p.MCP1.tDelay = 40;
p.MIP1a.kC = 0.010; p.MIP1a.tDelay = 30;  p.MIP1B.kC = 0.010; p.MIP1B.tDelay = 30;
p.plasmaBaseline = struct(); p.salivaBaseline = struct();
for nm = cytoNames
    p.plasmaBaseline.(char(nm)) = target.(char(nm)).baselinePlasma;
    p.salivaBaseline.(char(nm)) = target.(char(nm)).baseline;
end
end

function F = fit_unit(Tp, Ts, unit, C)
[target, cytoNames] = build_targets(Tp, Ts, unit, C.cytokineNames, C);
p = initial_params(target, cytoNames, C);
odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
optsFit = optimoptions('lsqnonlin', 'Display', 'off', 'MaxFunctionEvaluations', 6000, ...
    'MaxIterations', 1000, 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-10);

% ---- Stage A: shared upstream kinetics + TNF
q0A = [p.kLPS, p.kM, p.kM1, p.KLPS, p.nLPS, p.TNF.s, p.TNF.kmRNA, p.TNF.kC, p.TNF.kBS, p.TNF.kQ, p.TNF.kS];
lbA = [0.002 0.002 0.002 0.001 0.2   1e-4  0.01  0.0005  1e-6  0.001 0.001];
ubA = [0.5   0.5   0.5   5     6     500   3     0.3     0.5   0.5   0.5];
resA = @(q) stageA_residuals(q, [p.LPS0; 0; 0; 0; 0; 0; 0], odeOpts, target.TNF, p.plasmaBaseline.TNF, C);
[q, rn, ~, flag] = lsqnonlin(resA, q0A, lbA, ubA, optsFit);
p.kLPS = q(1); p.kM = q(2); p.kM1 = q(3); p.KLPS = q(4); p.nLPS = q(5);
p.TNF.s = q(6); p.TNF.kmRNA = q(7); p.TNF.kC = q(8); p.TNF.kBS = q(9); p.TNF.kQ = q(10); p.TNF.kS = q(11);
fprintf('  Stage A (upstream + TNF): resnorm = %.4g, exitflag = %d\n', rn, flag);

[tCore, YCore] = ode45(@(t, Y) coreODE(t, Y, p.kLPS, p.kM, p.kM1), C.tGrid, [p.LPS0; 0; 0], ...
    odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
LPSinterp = griddedInterpolant(tCore, YCore(:, 1), 'pchip', 'nearest');
M1interp  = griddedInterpolant(tCore, YCore(:, 3), 'pchip', 'nearest');

% ---- Stage B: every other cytokine separately
f = {'s', 'kmRNA', 'kC', 'tDelay', 'nDelay', 'kBS', 'kQ', 'kS', 'local', 'localDelay'};
for nm = cytoNames(cytoNames ~= "TNF")
    fn = char(nm);
    q0 = cellfun(@(x) p.(fn).(x), f);
    lb = [1e-4 0.01 0.0005   0 1  1e-6 0.001 0.001    0   0];
    ub = [500  3    0.3    300 10 0.5  0.5   0.5   3000 300];
    if ~C.hasLocalFun(nm), ub(9:10) = 0; end
    q0 = min(max(q0, lb), ub);
    resB = @(q) cytokine_residuals(q, odeOpts, p.KLPS, p.nLPS, LPSinterp, M1interp, ...
        target.(fn), p.plasmaBaseline.(fn), C.windowMap.(fn), C);
    [q, rn, ~, flag] = lsqnonlin(resB, q0, lb, ub, optsFit);
    for j = 1:numel(f), p.(fn).(f{j}) = q(j); end
    fprintf('  Stage B (%-6s): resnorm = %.4g, exitflag = %d\n', fn, rn, flag);
end

% ---- Stage C: all parameters jointly
[q0, lb, ub, parNames] = joint_vector(p, cytoNames, C.hasLocalFun, lbA, ubA);
x0 = zeros(3 + 4*numel(cytoNames), 1); x0(1) = p.LPS0;
tC = tic;
optsC = optimoptions(optsFit, 'MaxFunctionEvaluations', C.stageCMaxEval, 'MaxIterations', 3000, ...
    'UseParallel', C.useParallel, 'OutputFcn', @(x, ov, st) fit_progress(x, ov, st, tC));
resC = @(q) joint_residuals(q, p, x0, odeOpts, target, cytoNames, C);
[qJ, rn, residJ, flag, ~, ~, JJ] = lsqnonlin(resC, q0, lb, ub, optsC);
p = unpackJointParams(qJ, p, cytoNames, C.hasLocalFun);
fprintf('  Stage C (%d parameters): resnorm = %.4g, exitflag = %d (%.1f min)\n', numel(q0), rn, flag, toc(tC)/60);

F = struct('unit', unit, 'cytoNames', cytoNames, 'target', target, 'p', p, 'qJoint', qJ, ...
    'lb', lb, 'ub', ub, 'J', JJ, 'resid', residJ, 'parNames', parNames, 'resnorm', rn, 'exitflag', flag);
F.sim = simulate_model(p, cytoNames, C.tGrid, true);
F.metrics = derived_metrics(F.sim, p, cytoNames);
F.unc = [];
end

function [q0, lb, ub, names] = joint_vector(p, cytoNames, hasLocalFun, lbA, ubA)
q0 = [p.kLPS, p.kM, p.kM1, p.KLPS, p.nLPS];
lb = lbA(1:5); ub = ubA(1:5);
names = ["kLPS", "kM", "kM1", "KLPS", "nLPS"];
for nm = cytoNames
    c = p.(char(nm));
    isTNF = nm == "TNF";                            % TNF: no activation delay
    thisQ0 = [c.s, c.kmRNA, c.kC, c.tDelay, c.nDelay, c.kBS, c.kQ, c.kS];
    thisLB = [1e-4 0.01 0.0005 0 1 1e-6 0.001 0.001];
    thisUB = [500 3 0.3 300*~isTNF 1+9*~isTNF 0.5 0.5 0.5];
    thisNm = nm + "." + ["s", "kmRNA", "kC", "tDelay", "nDelay", "kBS", "kQ", "kS"];
    if hasLocalFun(nm)
        thisQ0 = [thisQ0, c.local, c.localDelay]; thisLB = [thisLB, 0, 0]; %#ok<AGROW>
        thisUB = [thisUB, 3000, 300]; thisNm = [thisNm, nm + ".local", nm + ".localDelay"]; %#ok<AGROW>
    end
    q0 = [q0, thisQ0]; lb = [lb, thisLB]; ub = [ub, thisUB]; names = [names, thisNm]; %#ok<AGROW>
end
q0 = min(max(q0, lb), ub);
end

function sim = simulate_model(p, cytoNames, tGrid, fine)
% Plasma, oral transit and saliva curves on tGrid (min).
nC = numel(cytoNames);
x0 = zeros(3 + 4*nC, 1); x0(1) = p.LPS0;
if fine, o = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 1.0);
else,    o = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 2.0); end
[tt, x] = ode45(@(tq, y) modelODE(tq, y, p, cytoNames), tGrid(:), x0, o);
if numel(tt) ~= numel(tGrid), error('ODE integration stopped early.'); end
sim.t = tt(:); sim.M1 = x(:, 3);
for k = 1:nC
    f = char(cytoNames(k)); i0 = 3 + (k-1)*4;
    sim.mrna(:, k)    = x(:, i0+1);
    sim.plasma(:, k)  = p.plasmaBaseline.(f) + x(:, i0+2);
    sim.transit(:, k) = x(:, i0+3);
    sim.saliva(:, k)  = p.salivaBaseline.(f) + x(:, i0+4);
end
end

function on = start_pool()
on = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
if on && isempty(gcp('nocreate')), parpool; end
end

function stop = fit_progress(~, ov, state, t0)
stop = false;
if strcmp(state, 'iter') && mod(ov.iteration, 10) == 0
    fprintf('    Stage C: iteration %4d, resnorm %.5g, %.1f min\n', ov.iteration, ov.resnorm, toc(t0)/60);
end
end


%% ========================================================================
%  DERIVED METRICS
%  ========================================================================
function names = metric_names()
% tDelay activation time; t12P ln2/kC; peakP/tPeakP, peakS/tPeakS peak and peak time;
% lag saliva minus plasma peak time; ratioSP saliva/plasma peak; s8overPeak saliva at
% 8 h / peak; lagC centroid lag 0-24 h; mtt mean transit time 1/kQ + 1/kS; sEnd saliva
% still rising at 24 h; foldP/foldS peak/baseline; rem8P/rem8S increase remaining at 8 h
% (% of the peak increase); halfP/halfS decline half-time; lagC8 centroid lag 0-8 h.
names = ["tDelay", "t12P", "peakP", "tPeakP", "peakS", "tPeakS", "lag", "ratioSP", "s8overPeak", ...
         "lagC", "mtt", "sEnd", "foldP", "foldS", "rem8P", "rem8S", "halfP", "halfS", "lagC8"];
end

function M = derived_metrics(sim, p, cytoNames)
% One row per cytokine, columns as metric_names; times in min.
nC = numel(cytoNames); M = nan(nC, numel(metric_names()));
t = sim.t(:); in8 = t <= 480;
for k = 1:nC
    f = char(cytoNames(k));
    yP = sim.plasma(:, k); yS = sim.saliva(:, k);
    bP = p.plasmaBaseline.(f); bS = p.salivaBaseline.(f);
    [pkP, tpP] = peak_time(t, yP);
    [pkS, tpS, atEnd] = peak_time(t, yS);
    yP8 = interp1(t, yP, 480); yS8 = interp1(t, yS, 480);
    M(k, :) = [p.(f).tDelay, log(2)/p.(f).kC, pkP, tpP, pkS, tpS, tpS - tpP, pkS/pkP, yS8/pkS, ...
        centroid(t, yS - bS) - centroid(t, yP - bP), 1/p.(f).kQ + 1/p.(f).kS, double(atEnd), ...
        pkP/bP, pkS/bS, 100*(yP8 - bP)/(pkP - bP), 100*(yS8 - bS)/(pkS - bS), ...
        decline_half(t, yP - bP, tpP), decline_half(t, yS - bS, tpS), ...
        centroid(t(in8), yS(in8) - bS) - centroid(t(in8), yP(in8) - bP)];
end
end

function [pk, tp, atEnd] = peak_time(t, y)
% Maximum, refined by a parabola through the highest grid point and its neighbours.
[pk, i] = max(y); tp = t(i); n = numel(y); atEnd = i == n;
if i > 1 && i < n
    y1 = y(i-1); y2 = y(i); y3 = y(i+1); den = y1 - 2*y2 + y3;
    if den < 0
        d = 0.5 * (y1 - y3) / den;
        tp = t(i) + d * (t(i+1) - t(i-1)) / 2;
        pk = y2 - 0.25 * (y1 - y3) * d;
    end
end
end

function h = decline_half(t, e, tPeak)
% Time from the peak until the increase e has fallen to half its maximum (NaN if never).
h = NaN; emax = max(e);
if ~(emax > 0), return; end
j = find(t > tPeak & e <= 0.5*emax, 1);
if isempty(j) || j == 1, return; end
t1 = t(j-1); t2 = t(j); e1 = e(j-1); e2 = e(j);
if e1 == e2, tc = t2; else, tc = t1 + (e1 - 0.5*emax) * (t2 - t1) / (e1 - e2); end
h = tc - tPeak;
end

function c = centroid(t, e)
e = max(e(:), 0); a = trapz(t, e);
if a <= 0, c = NaN; else, c = trapz(t, t .* e) / a; end
end


%% ========================================================================
%  BOOTSTRAP
%  ========================================================================
function U = uncertainty_bootstrap(Tp, Ts, unit, F, C, nBoot, maxIter, ciLevel, seed)
% Participants are resampled with replacement within each group, the targets rebuilt and
% stage C refitted from the point estimate.
rng(seed);
Pd = Tp(ismember(Tp.group, unit.groups), :);
Sd = Ts(ismember(Ts.group, unit.groups), :);
ids = unique([Pd.id; Sd.id]);
grpOfId = strings(numel(ids), 1); rowsP = cell(numel(ids), 1); rowsS = rowsP;
for i = 1:numel(ids)
    rowsP{i} = find(Pd.id == ids(i)); rowsS{i} = find(Sd.id == ids(i));
    if ~isempty(rowsP{i}), grpOfId(i) = Pd.group(rowsP{i}(1)); end
end
G = unique(grpOfId);
draws = zeros(nBoot, numel(ids));                   % drawn up front: independent of the workers
for b = 1:nBoot
    col = 0;
    for g = G(:)'
        ig = find(grpOfId == g); kk = numel(ig);
        draws(b, col+1:col+kk) = ig(randi(kk, kk, 1)); col = col + kk;
    end
end
cytoNames = F.cytoNames; nC = numel(cytoNames); nM = numel(metric_names());
tU = (0:2:C.t_end)'; nT = numel(tU);
Q = nan(nBoot, numel(F.qJoint)); SP = nan(nT, nC, nBoot); SS = SP; SPf = SP; SSf = SP;
MM = nan(nBoot, nC, nM); flags = nan(nBoot, 1);
p0 = F.p; q0 = F.qJoint; lb = F.lb; ub = F.ub;
nWorkers = 0; if C.useParallel, nWorkers = Inf; end
fprintf('  Bootstrap: %d resamples of %d participants...\n', nBoot, numel(ids));
parfor (b = 1:nBoot, nWorkers)
    try
        [Pb, Sb] = resample_tables(Pd, Sd, ids, rowsP, rowsS, draws(b, :));
        [tb, keep] = build_targets(Pb, Sb, unit, cytoNames, C);
        if numel(keep) < nC, error('a cytokine has no data in this resample'); end
        [qb, pb, fl] = refit_joint(p0, q0, lb, ub, tb, cytoNames, C, maxIter);
        flags(b) = fl;
        sb = simulate_model(pb, cytoNames, tU, false);
        Q(b, :) = qb; SP(:, :, b) = sb.plasma; SS(:, :, b) = sb.saliva;
        SPf(:, :, b) = sb.plasma ./ base_row(pb, cytoNames, 'plasma');   % fold change from the
        SSf(:, :, b) = sb.saliva ./ base_row(pb, cytoNames, 'saliva');   % resample's own baseline
        MM(b, :, :) = reshape(derived_metrics(sb, pb, cytoNames), [1 nC nM]);
    catch ME
        fprintf('    resample %d skipped: %s\n', b, ME.message);
    end
end
ok = all(isfinite(Q), 2);
nCap = sum(flags(ok) == 0);
fprintf('  Bootstrap: %d of %d resamples usable, %d stopped at the iteration limit.\n', sum(ok), nBoot, nCap);
a = 100*(1 - ciLevel)/2;
SP = SP(:, :, ok); SS = SS(:, :, ok); SPf = SPf(:, :, ok); SSf = SSf(:, :, ok); MM = MM(ok, :, :);
U = struct('t', tU, 'plasmaLo', pct(SP, a, 3), 'plasmaHi', pct(SP, 100 - a, 3), ...
    'salivaLo', pct(SS, a, 3), 'salivaHi', pct(SS, 100 - a, 3), ...
    'plasmaFoldLo', pct(SPf, a, 3), 'plasmaFoldHi', pct(SPf, 100 - a, 3), ...
    'salivaFoldLo', pct(SSf, a, 3), 'salivaFoldHi', pct(SSf, 100 - a, 3), 'metricDraws', MM, ...
    'metricLo', reshape(pct(MM, a, 1), nC, nM), 'metricHi', reshape(pct(MM, 100 - a, 1), nC, nM), ...
    'ciLevel', ciLevel, 'method', "bootstrap", 'nOK', sum(ok), 'n', nBoot, 'Q', Q(ok, :), 'nAtMaxIter', nCap);
end

function [Pb, Sb] = resample_tables(Pd, Sd, ids, rowsP, rowsS, ii)
% Rows of the drawn participants; a participant drawn twice gets two ids.
rp = cell(numel(ii), 1); rs = rp;
for k = 1:numel(ii)
    newId = ids(ii(k)) + "#" + k;
    A = Pd(rowsP{ii(k)}, :); A.id(:) = newId; rp{k} = A;
    B = Sd(rowsS{ii(k)}, :); B.id(:) = newId; rs{k} = B;
end
Pb = vertcat(rp{:}); Sb = vertcat(rs{:});
end

function [q, p, flag] = refit_joint(p, q0, lb, ub, target, cytoNames, C, maxIter)
% Stage C on resampled targets, started at the point estimate; baselines follow the data.
for nm = cytoNames
    p.salivaBaseline.(char(nm)) = target.(char(nm)).baseline;
    p.plasmaBaseline.(char(nm)) = target.(char(nm)).baselinePlasma;
end
x0 = zeros(3 + 4*numel(cytoNames), 1); x0(1) = p.LPS0;
opts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', maxIter*(numel(q0) + 1), 'FunctionTolerance', 1e-7, ...
    'StepTolerance', 1e-9, 'UseParallel', false);
rf = @(qq) joint_residuals(qq, p, x0, odeset('RelTol', 1e-6, 'AbsTol', 1e-8), target, cytoNames, C);
[q, ~, ~, flag] = lsqnonlin(rf, q0, lb, ub, opts);
p = unpackJointParams(q, p, cytoNames, C.hasLocalFun);
end

function b = base_row(p, cytoNames, which)
b = arrayfun(@(nm) p.([which 'Baseline']).(char(nm)), cytoNames(:)');
end

function Y = pct(X, p, dim)
% Percentile p (0-100) of X along dim, linear interpolation.
X = sort(X, dim); n = size(X, dim);
sz = size(X); sz(dim) = 1;
if n == 0, Y = nan(sz); return; end
pos = 1 + (n - 1) * p / 100; lo = floor(pos); hi = ceil(pos); w = pos - lo;
idx = repmat({':'}, 1, ndims(X));
idx{dim} = lo; A = X(idx{:});
idx{dim} = hi; B = X(idx{:});
Y = A + w .* (B - A);
end


%% ========================================================================
%  MODEL TABLES
%  ========================================================================
function T = model_tables(F, ciLevel)
mn = metric_names(); nC = numel(F.cytoNames); hasU = isfield(F, 'unc') && ~isempty(F.unc);

% all derived metrics with their confidence intervals (Table S2)
M = table(F.cytoNames(:), 'VariableNames', {'Cytokine'});
for j = 1:numel(mn)
    M.(char(mn(j))) = F.metrics(:, j);
    if hasU
        M.([char(mn(j)) '_lo']) = F.unc.metricLo(:, j); M.([char(mn(j)) '_hi']) = F.unc.metricHi(:, j);
    end
end

% peak and centroid lag, with a note where the lag is not resolved by the data
a = 100 * [(1 - ciLevel)/2, (1 + ciLevel)/2];
gm = @(m) F.metrics(:, mn == m);
note = strings(nC, 1);
for k = 1:nC
    j = find(F.parNames == F.cytoNames(k) + ".kQ" | F.parNames == F.cytoNames(k) + ".kS");
    atB = @(q) any(q(j) >= 0.99 * F.ub(j) & F.ub(j) > F.lb(j));
    if hasU && size(F.unc.metricDraws, 1) >= 10
        pEnd = mean(F.unc.metricDraws(:, k, mn == "sEnd"));
        pBound = mean(arrayfun(@(r) atB(F.unc.Q(r, :)), 1:size(F.unc.Q, 1)));
    else
        pEnd = F.metrics(k, mn == "sEnd"); pBound = double(atB(F.qJoint));
    end
    if pEnd > 0.2 || F.metrics(k, mn == "sEnd") == 1
        note(k) = "lower bound: saliva still rising at 24 h";
    elseif pBound > 0.2
        note(k) = "transfer at its fastest allowed rate: lag not resolved by the data";
    end
end
L = table(F.cytoNames(:), gm("lag"), gm("lagC"), note, 'VariableNames', {'Cytokine', 'lag', 'lagCentroid', 'note'});

% parameters (Tables S3-S4): estimate, bootstrap interval, allowed range, at a bound
lo = nan(numel(F.qJoint), 1); hi = lo;
if hasU, lo = pct(F.unc.Q, a(1), 1)'; hi = pct(F.unc.Q, a(2), 1)'; end
est = F.qJoint(:); lb = F.lb(:); ub = F.ub(:);
tol = 1e-9 * max(1, abs([lb ub]));
P = table(string(F.parNames(:)), est, lo, hi, lb, ub, double(est <= lb + tol(:, 1) | est >= ub - tol(:, 2)), ...
    'VariableNames', {'parameter', 'estimate', 'ci_lo', 'ci_hi', 'lower_bound', 'upper_bound', 'at_bound'});
cy = string(F.cytoNames(:));
B = table(cy, arrayfun(@(c) F.p.plasmaBaseline.(char(c)), cy), arrayfun(@(c) F.p.salivaBaseline.(char(c)), cy), ...
    'VariableNames', {'Cytokine', 'plasmaBaseline', 'salivaBaseline'});
T = struct('model_metrics', M, 'lags', L, 'parameters', P, 'baselines', B);
end
