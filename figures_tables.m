%% figures_tables.m
% Figures 2-6 and S2-S8 and the data and model tables of the manuscript, from the
% cleaned data (Tp, Ts, Tprot), the model settings C and the fit F made by
% endotoxemia_model.m. Run by run_all.m; FIGDIR and TABDIR are set there.
% The shift analysis bootstrap takes a few minutes.

opt = struct();
opt.groups   = F.unit.groups;                    % pooled groups
opt.allArms  = ["placebo", "anakinra", "tocilizumab", "dexamethason"];
opt.armLabel = ["Placebo", "Anakinra", "Tocilizumab", "Dexamethasone"];
opt.excl     = C.exclusionGroupMap;               % IL-1RA without anakinra, IL-6 without tocilizumab
opt.cyto     = C.cytokineNames;
opt.HL       = ["MIP1a", "MIP1B"];               % highlighted cytokines
opt.tMax     = 1440;                             % 0-24 h
opt.shifts   = 0:1:150;                          % pooled shift analysis, 1-min steps
opt.iShifts  = -60:1:150;                        % individual shift analysis
opt.minPtsI  = 4;                                % min paired points per participant
opt.nBootSh  = 2000;                             % bootstrap resamples for the shift CI
opt.nBootCI  = 2000;                             % bootstrap resamples for the CI of the data means
opt.seed     = 2026;
opt.ci       = 0.95;
opt.a        = 100 * [(1 - opt.ci)/2, (1 + opt.ci)/2];
opt.B        = struct('t1', 8.5, 't2', 23, 'gap', 0.9);     % broken time axis (h): 0-8.5 | 23-24
opt.xTick    = [0 2 4 6 8 24];
opt.xTickGrid = [0 4 8 24];
opt.cP = hs.C.ink; opt.cS = hs.C.accent;
opt.nCols = 5;
FULL = 174; MEDIUM = 120;                       % figure widths (mm)
for d = {FIGDIR, TABDIR}, if ~isfolder(d{1}), mkdir(d{1}); end, end
mn = metric_names(); nC = numel(F.cytoNames);
rng(opt.seed);
armCol = zeros(numel(opt.allArms), 3);
for a = 1:numel(opt.allArms), armCol(a, :) = arm_color(opt.allArms(a)); end

% data summaries (pooled, with the cytokine-specific exclusions), 0-24 h
SPp = cell(1, nC); SSp = SPp;
for k = 1:nC
    SPp{k} = summarise(Tp, opt.groups, F.cytoNames(k), opt.excl, 0, opt.tMax, opt.nBootCI, opt.ci);
    SSp{k} = summarise(Ts, opt.groups, F.cytoNames(k), opt.excl, 0, opt.tMax, opt.nBootCI, opt.ci);
end


%% ===== Figure 2: data per treatment group; arms_summary.csv ================
fig = hs.figure(FULL, 0.92);
tl = tiledlayout(fig, 4, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
dxA = [-6 -2 2 6];                                  % minutes, so error bars do not overlap
rawTab = table();
for k = 1:nC
    nm = opt.cyto(k);
    [iP, iS] = grid_index(k, opt.nCols);
    axP = nexttile(tl, iP); hold(axP, 'on');
    axS = nexttile(tl, iS); hold(axS, 'on');
    vP = []; vS = [];
    for a = 1:numel(opt.allArms)
        sp = summarise(Tp, opt.allArms(a), nm, [], -Inf, opt.tMax, opt.nBootCI, opt.ci);
        ss = summarise(Ts, opt.allArms(a), nm, [], -Inf, opt.tMax, opt.nBootCI, opt.ci);
        if isempty(sp.t) || isempty(ss.t), continue; end
        line_data(axP, sp, armCol(a, :), 's', dxA(a), opt.B);
        line_data(axS, ss, armCol(a, :), 'o', dxA(a), opt.B);
        vP = [vP; sp.lo; sp.hi]; vS = [vS; ss.lo; ss.hi]; %#ok<AGROW>
        rawTab = [rawTab; raw_row(sp, nm, "plasma", opt.allArms(a)); raw_row(ss, nm, "saliva", opt.allArms(a))]; %#ok<AGROW>
    end
    draw_lloq(axS, C, nm, opt.B);
    set(axP, 'YScale', 'log', 'YLim', log_lim(vP, 1.4));
    set(axS, 'YScale', 'log', 'YLim', log_lim(vS, 1.4));
    title(axP, cyto_label(nm));
    finish_pair(axP, axS, k, opt, -1.2);
end
xlabel(tl, 'Time after LPS (h)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
items = cell(numel(opt.allArms), 4);
for a = 1:numel(opt.allArms), items(a, :) = {char(opt.armLabel(a)), armCol(a, :), 'o', '-'}; end
grid_legend(tl, items);
hs.save(fig, fullfile(FIGDIR, 'Figure2'));
rawTab = add_vs_placebo(rawTab);
writetable(rawTab, fullfile(TABDIR, 'arms_summary.csv'));


%% ===== Shift analysis: Table 2, Tables S7 and S11, Figure S3 ================
grpOf = unique(Tp(:, {'id', 'group'}), 'rows');
T1 = table(); ROB = table(); curves = cell(1, nC);
for k = 1:nC
    nm = opt.cyto(k);
    [P, S] = pair_rows(Tp, Ts, opt.groups, nm, opt.excl, opt.tMax);
    [ids, ST] = shift_stats(P, S, opt.shifts);
    [r, nPairs] = r_from_stats(ST);
    nPart = sum(isfinite(ST(:, 1, 1)));
    [lo, hi] = fisher_ci(r, nPairs, nPart, opt.ci);
    [rOpt, jOpt] = max(r); sOpt = opt.shifts(jOpt);
    % bootstrap over participants, within arm, for the optimal shift
    g = strings(numel(ids), 1);
    for i = 1:numel(ids), g(i) = grpOf.group(find(grpOf.id == ids(i), 1)); end
    G = unique(g); bestB = nan(opt.nBootSh, 1); rOptB = bestB;
    for b = 1:opt.nBootSh
        idx = [];
        for u = G(:)'
            ig = find(g == u); idx = [idx; ig(randi(numel(ig), numel(ig), 1))]; %#ok<AGROW>
        end
        rb = r_from_stats(ST(idx, :, :));
        [rOptB(b), jb] = max(rb); bestB(b) = opt.shifts(jb);
    end
    ciS = round(pct(bestB, opt.a, 1));             % shifts are on a 1-min grid
    [loO, hiO] = fisher_ci(rOpt, nPairs(jOpt), nPart, opt.ci);
    curves{k} = struct('r', r, 'lo', lo, 'hi', hi, 'sOpt', sOpt, 'rOpt', rOpt, 'ciS', ciS);
    % within-participant Pearson r at no shift (individual fits, Figure 3)
    rInd = indiv_r(P, S, 0, 3);
    T1 = [T1; table(nm, nPart, nPairs(1), r(1), lo(1), hi(1), sOpt, ciS(1), ciS(2), rOpt, loO, hiO, ...
        median(rInd, 'omitnan'), quant(rInd, 0.25), quant(rInd, 0.75), ...
        'VariableNames', {'Cytokine', 'nParticipants', 'nPairs', 'r_noShift', 'r_noShift_lo', 'r_noShift_hi', ...
        'optimalShift_min', 'optimalShift_lo', 'optimalShift_hi', 'r_optimal', 'r_optimal_lo', 'r_optimal_hi', ...
        'medianIndividual_r', 'individual_r_q25', 'individual_r_q75'})]; %#ok<AGROW>

    % robustness: individual optimal shifts, without BLQ, protein-normalised
    best = nan(numel(ids), 1);
    for i = 1:numel(ids)
        rr = indiv_curve(P(P.id == ids(i), :), S(S.id == ids(i), :), opt.iShifts, opt.minPtsI);
        if all(isnan(rr)), continue; end
        [~, jb] = max(rr); best(i) = opt.iShifts(jb);
    end
    okB = ~isnan(best);
    nNeg = sum(best(okB) < 0); nZero = sum(best(okB) == 0); nPos = sum(best(okB) > 0);
    [~, STn] = shift_stats(P, S(~S.blq, :), opt.shifts); rN = r_from_stats(STn); [~, jN] = max(rN);
    Sp = add_protein(S, Tprot); Sp = Sp(isfinite(Sp.protein) & Sp.protein > 0, :);
    Sn = Sp; Sn.conc = Sp.conc ./ Sp.protein;
    [~, STr] = shift_stats(P, Sp, opt.shifts); rR = r_from_stats(STr);
    [~, STp] = shift_stats(P, Sn, opt.shifts); rPn = r_from_stats(STp); [~, jP] = max(rPn);
    [fRaw, ~] = geo_peak(Sp(Sp.time >= 0, :)); [fNrm, ~] = geo_peak(Sn(Sn.time >= 0, :));
    ROB = [ROB; table(nm, sum(okB), nNeg, nZero, nPos, sign_test(nPos, nNeg), median(best(okB)), ...
        quant(best(okB), 0.25), quant(best(okB), 0.75), opt.shifts(jN), sum(isfinite(STn(:, 1, 1))), ...
        fRaw, fNrm, rR(1), rPn(1), opt.shifts(jP), height(Sp), height(S), ...
        'VariableNames', {'Cytokine', 'nEvaluable', 'nNegative', 'nZero', 'nPositive', 'signTest_p', ...
        'medianIndividualShift', 'q25Shift', 'q75Shift', 'optimalShift_noBLQ', 'nPart_noBLQ', ...
        'foldRise_raw', 'foldRise_proteinNorm', 'r0_raw_withProtein', 'r0_proteinNorm', ...
        'optimalShift_proteinNorm', 'nSaliva_withProtein', 'nSaliva_all'})]; %#ok<AGROW>
end
writetable(T1, fullfile(TABDIR, 'shift_analysis.csv'));
writetable(ROB, fullfile(TABDIR, 'robustness.csv'));

% ---- Figure S3: r versus shift --------------------------------------------
fig = hs.figure(FULL, 0.48);
tl = tiledlayout(fig, 2, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nC
    cv = curves{k}; ax = nexttile(tl); hold(ax, 'on');
    yline(ax, 0, '-', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
    ok = isfinite(cv.r) & isfinite(cv.lo);
    if any(ok)
        patch(ax, cv.ciS([1 2 2 1]), [-0.5 -0.5 1 1], hs.C.accent, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');                     % 95% CI of the optimal shift
        hs.band(ax, opt.shifts(ok)', cv.lo(ok)', cv.hi(ok)', hs.C.ink);
        plot(ax, opt.shifts(ok), cv.r(ok), '-', 'Color', hs.C.ink, 'LineWidth', hs.LW_DATA);
        plot(ax, cv.sOpt, cv.rOpt, 'o', 'MarkerFaceColor', hs.C.ink, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
    end
    set(ax, 'XLim', [0 opt.shifts(end)], 'XTick', 0:50:150, 'YLim', [-0.5 1]);
    two_line_title(ax, cyto_label(opt.cyto(k)), sprintf('optimum %d min [%g%s%g]', ...
        cv.sOpt, cv.ciS(1), char(8211), cv.ciS(2)));
    hs.axes(ax);
    if mod(k - 1, opt.nCols) ~= 0, ax.YTickLabel = []; end
    if k <= opt.nCols, ax.XTickLabel = []; end
end
xlabel(tl, 'Shift of plasma relative to saliva (min)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
ylabel(tl, 'Within-participant correlation (r)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
hs.save(fig, fullfile(FIGDIR, 'FigureS3'));


%% ===== Figure 3 and Figure S2: matched saliva versus plasma ==================
for scMode = ["fits", "loops"]
    fig = hs.figure(FULL, 0.56);
    tl = tiledlayout(fig, 2, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    uT = unique(Ts.time(Ts.time <= opt.tMax)); tCol = hs.ordered(numel(uT));
    hKey = gobjects(numel(uT), 1);
    for k = 1:nC
        nm = opt.cyto(k);
        [P, S] = pair_rows(Tp, Ts, opt.groups, nm, opt.excl, opt.tMax);
        Q = pairs_at_shift(P, S, 0);                         % same pairs as Table 2
        ax = nexttile(tl); hold(ax, 'on');
        ids = unique(Q.id);
        if scMode == "fits"
            for i = 1:numel(ids)                            % per-participant fits
                Qi = Q(Q.id == ids(i), :);
                if height(Qi) < 3, continue; end
                cf = polyfit(Qi.lp, Qi.ls, 1); xx = [min(Qi.lp) max(Qi.lp)];
                plot(ax, 10.^xx, 10.^polyval(cf, xx), '-', 'Color', hs.tint(hs.C.context, 0.25), 'LineWidth', 0.5);
            end
        else
            for i = 1:numel(ids)                            % per-participant trajectories
                Qi = sortrows(Q(Q.id == ids(i), :), 'time');
                plot(ax, 10.^Qi.lp, 10.^Qi.ls, '-', 'Color', hs.tint(hs.C.context, 0.35), 'LineWidth', 0.4);
            end
        end
        for j = 1:numel(uT)
            m = Q.time == uT(j);
            if any(m)
                h = scatter(ax, 10.^Q.lp(m), 10.^Q.ls(m), 7, tCol(j, :), 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.2);
                if k == 1, hKey(j) = h; end
            end
        end
        if scMode == "fits"                                   % average (repeated-measures) fit
            [beta, alpha] = rm_line(Q);
            xx = [min(Q.lp) max(Q.lp)];
            plot(ax, 10.^xx, 10.^(alpha + beta*xx), '-', 'Color', hs.C.ink, 'LineWidth', 1.3);
            sub = sprintf('r = %.2f', T1.r_noShift(k));
        else                                                % mean loop
            [Gt, tg] = findgroups(Q.time);
            gx = splitapply(@mean, Q.lp, Gt); gy = splitapply(@mean, Q.ls, Gt);
            [~, o] = sort(tg);
            plot(ax, 10.^gx(o), 10.^gy(o), '-', 'Color', hs.C.ink, 'LineWidth', 1.1);
            plot(ax, 10.^gx(o(1)), 10.^gy(o(1)), 'o', 'MarkerFaceColor', hs.C.ink, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
            sub = sprintf('%d participants', numel(ids));
        end
        set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', robust_lim(10.^Q.lp), 'YLim', robust_lim(10.^Q.ls));
        two_line_title(ax, cyto_label(nm), sub);
        hs.axes(ax); decade_ticks(ax);
    end
    xlabel(tl, 'Plasma (pg/mL)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
    ylabel(tl, 'Saliva (pg/mL)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
    okKey = isgraphics(hKey);
    if any(okKey)
        keyLab = compose('%g h', round(uT(okKey)/60, 2));
        keyLab(uT(okKey) < 0) = {'pre-dose'};
        lgd = legend(hKey(okKey), keyLab, 'Box', 'off', 'FontSize', hs.FS_NOTE);
        lgd.Layout.Tile = 'east'; title(lgd, 'Time');
    end
    if scMode == "fits"
        hs.save(fig, fullfile(FIGDIR, 'Figure3'));
    else
        hs.save(fig, fullfile(FIGDIR, 'FigureS2'));
    end
end


%% ===== Figure 4: fold change from baseline =================================
fig = hs.figure(FULL, 0.56);
tl = tiledlayout(fig, 2, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nC
    nm = F.cytoNames(k); f = char(nm);
    bP = F.p.plasmaBaseline.(f); bS = F.p.salivaBaseline.(f);
    ax = nexttile(tl); hold(ax, 'on');
    yline(ax, 1, '-', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
    model_band(ax, F, k, "plasma", "fold", opt.cP, opt.B); model_band(ax, F, k, "saliva", "fold", opt.cS, opt.B);
    model_line(ax, F, k, "plasma", "fold", opt.cP, '-', opt.B); model_line(ax, F, k, "saliva", "fold", opt.cS, '-', opt.B);
    point_data(ax, SPp{k}, opt.cP, 's', bP, opt.B); point_data(ax, SSp{k}, opt.cS, 'o', bS, opt.B);
    vals = [SPp{k}.lo/bP; SPp{k}.hi/bP; SSp{k}.lo/bS; SSp{k}.hi/bS; F.sim.plasma(:, k)/bP; F.sim.saliva(:, k)/bS];
    set(ax, 'YScale', 'log', 'YLim', log_lim(vals, 1.3));
    two_line_title(ax, cyto_label(nm), sprintf('P %s | S %s', fmt_conc(bP), fmt_conc(bS)));
    hs.axes(ax); bx_axis(ax, opt, 0, opt.xTickGrid); bx_break(ax, opt.B);
    if k <= opt.nCols, ax.XTickLabel = []; end
end
xlabel(tl, 'Time after LPS (h)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
ylabel(tl, 'Fold change from baseline', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
grid_legend(tl, {'plasma', opt.cP, 's', '-'; 'saliva', opt.cS, 'o', '-'});
hs.save(fig, fullfile(FIGDIR, 'Figure4'));


%% ===== Figure 5: model curves scaled to their own peak ======================
fig = hs.figure(FULL, 0.48);
tl = tiledlayout(fig, 2, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
tH = F.sim.t/60;
for k = 1:nC
    nm = F.cytoNames(k); f = char(nm);
    ax = nexttile(tl); hold(ax, 'on');
    yP = F.sim.plasma(:, k) - F.p.plasmaBaseline.(f);
    yS = F.sim.saliva(:, k) - F.p.salivaBaseline.(f);
    yQ = F.sim.transit(:, k);
    plot(ax, bx(tH, opt.B), yQ / max(max(yQ), eps), '--', 'Color', hs.C.context, 'LineWidth', hs.LW_THIN + 0.2);
    plot(ax, bx(tH, opt.B), yP / max(max(yP), eps), '-', 'Color', opt.cP, 'LineWidth', hs.LW_DATA);
    plot(ax, bx(tH, opt.B), yS / max(max(yS), eps), '-', 'Color', opt.cS, 'LineWidth', hs.LW_DATA);
    two_line_title(ax, cyto_label(nm), sprintf('lag %.0f min', F.metrics(k, mn == "lag")));
    set(ax, 'YLim', [0 1.08], 'YTick', [0 0.5 1]);
    hs.axes(ax); bx_axis(ax, opt, 0); bx_break(ax, opt.B);
    if mod(k - 1, opt.nCols) ~= 0, ax.YTickLabel = []; end
    if k <= opt.nCols, ax.XTickLabel = []; end
end
xlabel(tl, 'Time after LPS (h)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
ylabel(tl, 'Response, fraction of own peak', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
grid_legend(tl, {'plasma', opt.cP, 'none', '-'; 'oral transit', hs.C.context, 'none', '--'; 'saliva', opt.cS, 'none', '-'});
hs.save(fig, fullfile(FIGDIR, 'Figure5'));


%% ===== Table 3, Table S9 and Figure 6: kinetics =============================
D = F.unc.metricDraws;                               % nDraw x nC x nMetric
ciM = @(k, name) pct(D(:, k, mn == name), opt.a, 1);
T2 = table(); TH = table();
for k = 1:nC
    nm = F.cytoNames(k);
    g = @(name) F.metrics(k, mn == name);
    c = @(name) reshape(ciM(k, name), 1, 2);
    T2 = [T2; table(nm, g("tPeakP"), g("tPeakS"), g("lag"), c("lag"), g("foldP"), c("foldP"), ...
        g("foldS"), c("foldS"), g("rem8P"), c("rem8P"), g("rem8S"), c("rem8S"), ...
        'VariableNames', {'Cytokine', 'plasmaPeak_min', 'salivaPeak_min', 'lag_min', 'lag_CI', ...
        'plasmaFoldRise', 'plasmaFoldRise_CI', 'salivaFoldRise', 'salivaFoldRise_CI', ...
        'plasmaRemaining8h_pct', 'plasmaRemaining8h_CI', 'salivaRemaining8h_pct', 'salivaRemaining8h_CI'})]; %#ok<AGROW>
    TH = [TH; table(nm, g("t12P"), c("t12P"), g("halfP"), c("halfP"), g("halfS"), c("halfS"), ...
        g("mtt"), c("mtt"), g("lagC"), c("lagC"), g("lagC8"), c("lagC8"), ...
        'VariableNames', {'Cytokine', 'plasma_t12_ln2_over_kC', 'plasma_t12_CI', 'plasmaDeclineHalfTime', ...
        'plasmaDeclineHalfTime_CI', 'salivaDeclineHalfTime', 'salivaDeclineHalfTime_CI', 'meanTransitTime', ...
        'meanTransitTime_CI', 'centroidLag_0_24h', 'centroidLag_0_24h_CI', 'centroidLag_0_8h', 'centroidLag_0_8h_CI'})]; %#ok<AGROW>
end
T2 = sortrows(T2, 'lag_min');
writetable(T2, fullfile(TABDIR, 'kinetics.csv'));
writetable(TH, fullfile(TABDIR, 'half_lives.csv'));

% Figure 6: A plasma t1/2, B saliva t1/2, C lag, D fold-rise ratio, E remaining ratio
panels = {{'Plasma t_{1/2}', '(min)'}, @(M) M(:, mn == "halfP"), true, [30 50 100 200 500];
          {'Saliva t_{1/2}', '(min)'}, @(M) M(:, mn == "halfS"), true, [30 50 100 200 500 1000];
          {'Saliva lag', '(min)'},     @(M) M(:, mn == "lag"),   false, 0:100:300;
          {'Fold rise', 'saliva / plasma'}, @(M) M(:, mn == "foldS") ./ M(:, mn == "foldP"), true, [0.001 0.01 0.1 1];
          {'Rise left at 8 h', '(% of peak rise)'}, [], true, [0.1 1 10 100]};   % E: plasma and saliva side by side
[~, ord] = sort(F.metrics(:, mn == "lag")); rowCyto = F.cytoNames(ord);
fig = hs.figure(FULL, 0.46);
tl = tiledlayout(fig, 1, size(panels, 1), 'TileSpacing', 'compact', 'Padding', 'compact');
ax6 = gobjects(size(panels, 1), 1);
for j = 1:size(panels, 1)
    ax = nexttile(tl); hold(ax, 'on'); ax6(j) = ax;
    for r = 1:numel(rowCyto)
        k = find(F.cytoNames == rowCyto(r));
        cc = opt.cP; if ismember(rowCyto(r), opt.HL), cc = opt.cS; end
        if isempty(panels{j, 2})            % panel E: plasma (square) and saliva (circle), joined
            vv = [F.metrics(k, mn == "rem8P"), F.metrics(k, mn == "rem8S")];
            vv = max(vv, 0.01);             % floor for the log axis (values below 0.01% shown at 0.01)
            dy = [-0.13 0.13];
            plot(ax, vv, r + dy, '-', 'Color', hs.tint(cc, 0.45), 'LineWidth', 0.8);
            for q = 1:2
                nmq = ["rem8P", "rem8S"]; mk = 'so';
                ci = max(pct(D(:, k, mn == nmq(q)), opt.a, 1), 0.01);
                plot(ax, ci, [r r] + dy(q), '-', 'Color', cc, 'LineWidth', 0.8);
                plot(ax, vv(q), r + dy(q), mk(q), 'MarkerFaceColor', cc, 'MarkerEdgeColor', 'w', 'MarkerSize', 4.5);
            end
            continue
        end
        v = panels{j, 2}(F.metrics(k, :));
        dr = panels{j, 2}(reshape(D(:, k, :), [], numel(mn)));
        ci = pct(dr(isfinite(dr)), opt.a, 1);
        plot(ax, ci, [r r], '-', 'Color', cc, 'LineWidth', 1);
        plot(ax, v, r, 'o', 'MarkerFaceColor', cc, 'MarkerEdgeColor', 'w', 'MarkerSize', 4.5);
    end
    if panels{j, 3}, set(ax, 'XScale', 'log'); end
    xl = xlim(ax); tk = panels{j, 4}; tk = tk(tk >= xl(1) & tk <= xl(2));
    if numel(tk) >= 2, set(ax, 'XTick', tk, 'XTickLabel', compose('%g', tk), 'XTickLabelRotation', 0); end
    if j == 4, xline(ax, 1, '-', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN); end
    if j == 5
        xlim(ax, [0.01 100]); set(ax, 'XTick', [0.01 0.1 1 10 100], 'XTickLabel', {'<0.01', '0.1', '1', '10', '100'});
        hp = plot(ax, NaN, NaN, 's', 'MarkerFaceColor', opt.cP, 'MarkerEdgeColor', 'w', 'MarkerSize', 4.5);
        hsl = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', opt.cP, 'MarkerEdgeColor', 'w', 'MarkerSize', 4.5);
        legend(ax, [hp hsl], {'plasma', 'saliva'}, 'Box', 'off', 'Location', 'northwest', 'FontSize', hs.FS_NOTE);
    end
    set(ax, 'YTick', 1:numel(rowCyto), 'YDir', 'reverse', 'YLim', [0.4 numel(rowCyto) + 0.6]);
    if j == 1
        set(ax, 'YTickLabel', arrayfun(@(nm) cyto_label(nm), rowCyto, 'UniformOutput', false));
    else
        set(ax, 'YTickLabel', []);
    end
    xlabel(ax, panels{j, 1});
    hs.axes(ax);
end
hs.panelLabels(ax6);
hs.save(fig, fullfile(FIGDIR, 'Figure6'));


%% ===== Table S6: spread of data and model, plasma versus saliva =============
SPR = table();
iT = @(U, t) find(abs(U - t) == min(abs(U - t)), 1);
for k = 1:nC
    nm = F.cytoNames(k); row = table(nm, 'VariableNames', {'Cytokine'});
    for mat = ["plasma", "saliva"]
        if mat == "plasma", Sd = SPp{k}; else, Sd = SSp{k}; end
        post = Sd.t >= 0;
        gsd = Sd.gsd(post); dW = Sd.hi(post) ./ Sd.lo(post);
        [~, iPk] = max(Sd.m(post)); tt = Sd.t(post);
        lo = F.unc.([char(mat) 'Lo'])(:, k); hi = F.unc.([char(mat) 'Hi'])(:, k);
        mW = hi ./ lo;
        mWd = arrayfun(@(t) mW(iT(F.unc.t, t)), tt);           % model width at the data times
        tpk = F.metrics(k, mn == "tPeak" + upper(extractBefore(mat, 2)));
        pre = char(mat);
        row.([pre '_GSD_median'])          = median(gsd, 'omitnan');
        row.([pre '_GSD_atPeak'])          = gsd(iPk);
        row.([pre '_dataCI_foldWidth_med']) = median(dW, 'omitnan');
        row.([pre '_modelCI_foldWidth_med']) = median(mWd, 'omitnan');
        row.([pre '_modelCI_foldWidth_atPeak']) = mW(iT(F.unc.t, tpk));
    end
    row.saliva_over_plasma_GSD = row.saliva_GSD_median / row.plasma_GSD_median;
    row.saliva_over_plasma_modelCI = row.saliva_modelCI_foldWidth_med / row.plasma_modelCI_foldWidth_med;
    SPR = [SPR; row]; %#ok<AGROW>
end
writetable(SPR, fullfile(TABDIR, 'spread.csv'));


%% ===== Table S1 (saliva), Table S8 and Figure S7: fit quality ===============
BLQ = table(); tS = unique(Ts.time(Ts.time >= 0 & Ts.time <= opt.tMax));
for k = 1:nC
    R = data_rows(Ts, opt.groups, opt.cyto(k), opt.excl); R = R(R.time >= 0 & R.time <= opt.tMax, :);
    row = table(opt.cyto(k), height(R), 100*mean(R.blq), 'VariableNames', {'Cytokine', 'n', 'pctBLQ_all'});
    for t = tS'
        m = R.time == t; v = NaN; if any(m), v = 100*mean(R.blq(m)); end
        row.(sprintf('t%d_min', round(t))) = v;
    end
    BLQ = [BLQ; row]; %#ok<AGROW>
end
writetable(BLQ, fullfile(TABDIR, 'saliva_below_lloq.csv'));

PK = table();
for k = 1:nC
    d = F.target.(char(F.cytoNames(k)));
    [dpP, iP] = max(d.mPlasma); [dpS, iS] = max(d.mSaliva);
    PK = [PK; table(F.cytoNames(k), F.metrics(k, mn == "peakP"), dpP, 100*F.metrics(k, mn == "peakP")/dpP, ...
        F.metrics(k, mn == "tPeakP") - d.tPlasma(iP), F.metrics(k, mn == "peakS"), dpS, ...
        100*F.metrics(k, mn == "peakS")/dpS, F.metrics(k, mn == "tPeakS") - d.tSaliva(iS), ...
        'VariableNames', {'Cytokine', 'plasmaModel', 'plasmaObserved', 'plasmaPct', 'plasmaDeltaPeakTime', ...
        'salivaModel', 'salivaObserved', 'salivaPct', 'salivaDeltaPeakTime'})]; %#ok<AGROW>
end
writetable(PK, fullfile(TABDIR, 'peaks.csv'));

fig = hs.figure(MEDIUM, 0.52);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axO = gobjects(2, 1); mats = ["plasma", "saliva"]; GM = zeros(2, 2);
for j = 1:2
    ax = nexttile(tl); hold(ax, 'on'); axO(j) = ax;
    allO = []; allM = [];
    order = [find(~ismember(F.cytoNames, opt.HL)), find(ismember(F.cytoNames, opt.HL))];
    for k = order
        d = F.target.(char(F.cytoNames(k)));
        if j == 1, tt = d.tPlasma; ob = d.mPlasma; y = F.sim.plasma(:, k);
        else,      tt = d.tSaliva; ob = d.mSaliva; y = F.sim.saliva(:, k); end
        pr = interp1(F.sim.t, y, tt, 'pchip');
        cc = hs.C.ink; if ismember(F.cytoNames(k), opt.HL), cc = hs.C.accent; end
        plot(ax, pr, ob, 'o', 'MarkerFaceColor', cc, 'MarkerEdgeColor', 'w', 'MarkerSize', 3.5);
        allO = [allO; ob(:)]; allM = [allM; pr(:)]; %#ok<AGROW>
    end
    lim = log_lim([allO; allM], 1.6);
    plot(ax, lim, lim, '-', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
    plot(ax, lim, 2*lim, '--', lim, lim/2, '--', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
    e = abs(log(allO ./ allM)); GM(j, :) = [exp(mean(e)), 100*mean(e <= log(2))];
    corner_note(ax, sprintf('GMFE %.2f\n%.0f%% within 2-fold', GM(j, 1), GM(j, 2)), 'nw', hs.C.ink2);
    if j == 1, corner_note(ax, ['MIP-1' char(945) ', MIP-1' char(946)], 'se', hs.C.accent); end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lim, 'YLim', lim); axis(ax, 'square'); decade_ticks(ax);
    xlabel(ax, 'Fitted (pg/mL)'); ylabel(ax, 'Observed mean (pg/mL)');
    title(ax, [upper(extractBefore(char(mats(j)), 2)) extractAfter(char(mats(j)), 1)]);
    hs.axes(ax);
end
hs.panelLabels(axO);
hs.save(fig, fullfile(FIGDIR, 'FigureS7'));


%% ===== Figure S6 (log axes) and Figure S5 (linear axes) =====================
for sc = ["log", "linear"]
    fig = hs.figure(FULL, 0.88);
    tl = tiledlayout(fig, 4, opt.nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:nC
        nm = F.cytoNames(k);
        [iP, iS] = grid_index(k, opt.nCols);
        axP = nexttile(tl, iP); hold(axP, 'on');
        axS = nexttile(tl, iS); hold(axS, 'on');
        model_band(axP, F, k, "plasma", "abs", opt.cP, opt.B); model_band(axS, F, k, "saliva", "abs", opt.cS, opt.B);
        if sc == "linear"
            draw_subjects(axP, SPp{k}, opt.cP, opt.B); draw_subjects(axS, SSp{k}, opt.cS, opt.B);
        end
        model_line(axP, F, k, "plasma", "abs", opt.cP, '-', opt.B); model_line(axS, F, k, "saliva", "abs", opt.cS, '-', opt.B);
        draw_lloq(axS, C, nm, opt.B);
        point_data(axP, SPp{k}, opt.cP, 's', 1, opt.B); point_data(axS, SSp{k}, opt.cS, 'o', 1, opt.B);
        if sc == "log"
            set(axP, 'YScale', 'log', 'YLim', log_lim([SPp{k}.lo; SPp{k}.hi; F.sim.plasma(:, k)], 1.5));
            set(axS, 'YScale', 'log', 'YLim', log_lim([SSp{k}.lo; SSp{k}.hi; F.sim.saliva(:, k)], 1.5));
        else
            ylim(axP, [0 max(ylim(axP))]);
        end
        title(axP, cyto_label(nm));
        finish_pair(axP, axS, k, opt, 0);
    end
    xlabel(tl, 'Time after LPS (h)', 'FontSize', hs.FS_LABEL, 'FontName', hs.FONT);
    grid_legend(tl, {'plasma', opt.cP, 's', '-'; 'saliva', opt.cS, 'o', '-'});
    if sc == "log"
        hs.save(fig, fullfile(FIGDIR, 'FigureS6'));
    else
        hs.save(fig, fullfile(FIGDIR, 'FigureS5'));
    end
end


%% ===== Figure S4: worked example (MIP-1alpha) ===============================
k = find(F.cytoNames == "MIP1a"); f = 'MIP1a';
if ~isempty(k)
    fig = hs.figure(MEDIUM, 0.80);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl); hold(ax1, 'on');
    model_band(ax1, F, k, "plasma", "abs", opt.cP, opt.B); draw_subjects(ax1, SPp{k}, opt.cP, opt.B);
    model_line(ax1, F, k, "plasma", "abs", opt.cP, '-', opt.B); point_data(ax1, SPp{k}, opt.cP, 's', 1, opt.B);
    ylabel(ax1, 'Plasma (pg/mL)'); title(ax1, ['MIP-1' char(945) ' plasma']);
    hs.axes(ax1); bx_axis(ax1, opt, 0); ax1.XTickLabel = []; bx_break(ax1, opt.B);
    ax2 = nexttile(tl); hold(ax2, 'on');
    model_band(ax2, F, k, "saliva", "abs", opt.cS, opt.B); draw_subjects(ax2, SSp{k}, opt.cS, opt.B);
    model_line(ax2, F, k, "saliva", "abs", opt.cS, '-', opt.B);
    plot(ax2, bx(F.sim.t/60, opt.B), F.p.salivaBaseline.(f) + F.sim.transit(:, k), '--', 'Color', hs.C.ink2, 'LineWidth', hs.LW_THIN);
    yline(ax2, F.p.salivaBaseline.(f), ':', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
    point_data(ax2, SSp{k}, opt.cS, 'o', 1, opt.B);
    xlabel(ax2, 'Time after LPS (h)'); ylabel(ax2, 'Saliva (pg/mL)'); title(ax2, ['MIP-1' char(945) ' saliva']);
    hs.axes(ax2); bx_axis(ax2, opt, 0); bx_break(ax2, opt.B);
    hs.panelLabels([ax1 ax2]);
    hs.save(fig, fullfile(FIGDIR, 'FigureS4'));
end


%% ===== Figure S8: salivary total protein, 0-24 h =============================
fig = hs.figure(MEDIUM, 0.55);
ax = axes(fig); hold(ax, 'on');
sets = {opt.groups, "dexamethason"}; labs = {'Placebo + anakinra + tocilizumab', 'Dexamethasone'};
cols = [hs.C.ink; arm_color("dexamethason")]; hP = gobjects(2, 1);
for jj = 1:2
    gid = unique(Ts.id(ismember(Ts.group, sets{jj})));
    Rp = Tprot(ismember(Tprot.id, gid) & Tprot.time >= 0 & Tprot.time <= opt.tMax & Tprot.protein > 0, :);
    if isempty(Rp), continue; end
    [Gp, tvp] = findgroups(Rp.time);
    m = splitapply(@(v) mean(log(v)), Rp.protein, Gp);
    se = splitapply(@(v) std(log(v)) / sqrt(numel(v)), Rp.protein, Gp);
    ref = m(1); dx = (jj - 1.5) * 4 / 60;
    xx = bx(tvp/60 + dx, opt.B);
    errorbar(ax, xx, exp(m - ref), exp(m - ref) - exp(m - ref - 1.96*se), exp(m - ref + 1.96*se) - exp(m - ref), ...
        'o', 'Color', cols(jj, :), 'MarkerFaceColor', cols(jj, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 4.5, ...
        'CapSize', 0, 'LineWidth', 0.8);
    hP(jj) = plot(ax, xx, exp(m - ref), '-', 'Color', cols(jj, :), 'LineWidth', hs.LW_DATA);
end
yline(ax, 1, ':', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN);
xlabel(ax, 'Time after LPS (h)'); ylabel(ax, 'Salivary total protein (fold of 0.5 h)');
legend(ax, hP(isgraphics(hP)), labs(isgraphics(hP)), 'Box', 'off', 'Location', 'northeast', 'FontSize', hs.FS_NOTE);
hs.axes(ax); bx_axis(ax, opt, 0); bx_break(ax, opt.B);
hs.save(fig, fullfile(FIGDIR, 'FigureS8'));


%% =============================================================================
%  LOCAL FUNCTIONS
%  =============================================================================
function names = metric_names()
    names = ["tDelay", "t12P", "peakP", "tPeakP", "peakS", "tPeakS", "lag", "ratioSP", "s8overPeak", ...
             "lagC", "mtt", "sEnd", "foldP", "foldS", "rem8P", "rem8S", "halfP", "halfS", "lagC8"];
end

% ---- data selection -----------------------------------------------------------
function R = data_rows(T, groups, cyto, excl)
    % Day-0 rows of one cytokine for the given arms, with the cytokine-specific
    % exclusion (excl: containers.Map cytokine -> arm, or [] for none).
    R = T(ismember(T.group, string(groups)) & T.cyto == cyto, :);
    key = char(cyto);
    if ~isempty(excl) && isKey(excl, key)
        R = R(R.group ~= string(excl(key)), :);
    end
end

function [P, S] = pair_rows(Tp, Ts, groups, cyto, excl, tMax)
    % Plasma and saliva rows for the shift analysis: pre-dose to tMax.
    P = data_rows(Tp, groups, cyto, excl); P = P(P.conc > 0 & P.time <= tMax, :);
    S = data_rows(Ts, groups, cyto, excl); S = S(S.conc > 0 & S.time <= tMax, :);
end

function S = summarise(T, groups, cyto, excl, tMin, tMax, nBoot, ciLevel)
    % Geometric mean per nominal time with a percentile-bootstrap CI over
    % participants, the between-participant geometric SD, and each series.
    R = data_rows(T, groups, cyto, excl);
    R = R(R.time >= tMin & R.time <= tMax & R.conc > 0, :);
    S = struct('t', [], 'm', [], 'lo', [], 'hi', [], 'n', [], 'gsd', [], 'subj', {{}});
    if isempty(R), return; end
    [G, tv] = findgroups(R.time);
    nT = numel(tv); S.t = tv; S.m = nan(nT, 1); S.lo = S.m; S.hi = S.m; S.n = S.m; S.gsd = S.m;
    a = (1 - ciLevel)/2;
    for i = 1:nT
        v = R.conc(G == i); n = numel(v); lv = log(v);
        S.n(i) = n; S.m(i) = exp(mean(lv));
        if n >= 2
            S.gsd(i) = exp(std(lv));
            bm = sort(exp(mean(lv(randi(n, n, nBoot)), 1)));
            S.lo(i) = bm(max(1, round(a*nBoot))); S.hi(i) = bm(min(nBoot, round((1 - a)*nBoot)));
        else
            S.lo(i) = S.m(i); S.hi(i) = S.m(i);
        end
    end
    ids = unique(R.id); S.subj = cell(numel(ids), 1);
    for k = 1:numel(ids)
        Rk = sortrows(R(R.id == ids(k), :), 'time'); S.subj{k} = [Rk.time, Rk.conc];
    end
end

function row = raw_row(S, cyto, mat, arm)
    % Per arm: pre-dose, 0.5-h value, peak (geometric mean) and its time, fold rise vs 0.5 h, 24-h value.
    pre = S.m(S.t < 0); if isempty(pre), pre = NaN; else, pre = pre(end); end
    i05 = find(S.t == 30, 1); b = NaN; if ~isempty(i05), b = S.m(i05); end
    post = S.t >= 0; tt = S.t(post); mm = S.m(post);
    [pk, ip] = max(mm);
    i24 = find(S.t == 1440, 1); v24 = NaN; if ~isempty(i24), v24 = S.m(i24); end
    row = table(string(cyto), string(mat), string(arm), pre, b, pk, tt(ip), pk / b, v24 / b, ...
        'VariableNames', {'Cytokine', 'matrix', 'arm', 'preDose', 'at30min', 'peak', 'peakTime_min', ...
        'foldRise_vs30min', 'at24h_vs30min'});
end

function T = add_vs_placebo(T)
    % Peak and fold rise of each arm relative to placebo (same cytokine and matrix).
    T.peak_vs_placebo = nan(height(T), 1); T.foldRise_vs_placebo = nan(height(T), 1);
    for i = 1:height(T)
        j = find(T.Cytokine == T.Cytokine(i) & T.matrix == T.matrix(i) & T.arm == "placebo", 1);
        if ~isempty(j)
            T.peak_vs_placebo(i) = T.peak(i) / T.peak(j);
            T.foldRise_vs_placebo(i) = T.foldRise_vs30min(i) / T.foldRise_vs30min(j);
        end
    end
end

function S = add_protein(S, Tprot)
    % Salivary total protein of the same sample (NaN if missing).
    kS = S.id + "|" + string(round(S.time));
    kP = Tprot.id + "|" + string(round(Tprot.time));
    [has, loc] = ismember(kS, kP);
    S.protein = nan(height(S), 1); S.protein(has) = Tprot.protein(loc(has));
end

function [fold, tPeak] = geo_peak(R)
    fold = NaN; tPeak = NaN;
    if isempty(R), return; end
    [G, tv] = findgroups(R.time);
    gm = splitapply(@(v) exp(mean(log(v))), R.conc, G);
    [mx, i] = max(gm); fold = mx / gm(1); tPeak = tv(i);
end

% ---- shift analysis -----------------------------------------------------------
function [ids, ST] = shift_stats(P, S, shifts)
    % Per participant and shift: centred sums Sxx, Syy, Sxy and n of
    % log10 saliva(t) and log10 plasma(t - shift), plasma linearly
    % interpolated within the participant. The repeated-measures correlation
    % (Bakdash & Marusich 2017) at a shift is sum(Sxy)/sqrt(sum(Sxx)sum(Syy)),
    % so resampling participants only needs these sums.
    ids = intersect(unique(P.id), unique(S.id));
    ST = nan(numel(ids), numel(shifts), 4);
    for i = 1:numel(ids)
        Pi = P(P.id == ids(i), :); Si = S(S.id == ids(i), :);
        [tp, ia] = unique(Pi.time);
        if numel(tp) < 2, continue; end
        yp = log10(Pi.conc(ia));
        for j = 1:numel(shifts)
            tq = max(Si.time - shifts(j), tp(1));       % before the first sample: pre-dose value
            ok = tq <= tp(end);
            if sum(ok) < 3, continue; end
            x = log10(Si.conc(ok)); y = interp1(tp, yp, tq(ok));
            x = x - mean(x); y = y - mean(y);
            ST(i, j, :) = [sum(x.*x), sum(y.*y), sum(x.*y), sum(ok)];
        end
    end
end

function [r, n] = r_from_stats(ST)
    sxx = sum(ST(:, :, 1), 1, 'omitnan'); syy = sum(ST(:, :, 2), 1, 'omitnan');
    sxy = sum(ST(:, :, 3), 1, 'omitnan'); n = sum(ST(:, :, 4), 1, 'omitnan');
    r = sxy ./ sqrt(sxx .* syy);
    r(~isfinite(r)) = NaN;
end

function [lo, hi] = fisher_ci(r, n, nSubj, ciLevel)
    dfr = n - nSubj - 1;
    zc = sqrt(2) * erfinv(ciLevel);
    se = 1 ./ sqrt(max(dfr - 1, 1));
    z = atanh(max(min(r, 0.999999), -0.999999));
    lo = tanh(z - zc*se); hi = tanh(z + zc*se);
end

function Q = pairs_at_shift(P, S, shift)
    % Matched pairs (log10 plasma interpolated at the salivary time - shift).
    ids = intersect(unique(P.id), unique(S.id));
    id = strings(0, 1); tm = []; lp = []; ls = [];
    for i = 1:numel(ids)
        Pi = P(P.id == ids(i), :); Si = S(S.id == ids(i), :);
        [tp, ia] = unique(Pi.time);
        if numel(tp) < 2, continue; end
        tq = max(Si.time - shift, tp(1)); ok = tq <= tp(end);
        if sum(ok) < 3, continue; end
        id = [id; repmat(ids(i), sum(ok), 1)]; tm = [tm; Si.time(ok)]; %#ok<AGROW>
        lp = [lp; interp1(tp, log10(Pi.conc(ia)), tq(ok))]; ls = [ls; log10(Si.conc(ok))]; %#ok<AGROW>
    end
    Q = table(id, tm, lp, ls, 'VariableNames', {'id', 'time', 'lp', 'ls'});
end

function [beta, alpha] = rm_line(Q)
    % Common within-participant slope (saliva on plasma) and mean intercept.
    ids = unique(Q.id); sxy = 0; sxx = 0; mx = nan(numel(ids), 1); my = mx;
    for i = 1:numel(ids)
        m = Q.id == ids(i); x = Q.lp(m); y = Q.ls(m);
        mx(i) = mean(x); my(i) = mean(y);
        sxy = sxy + sum((x - mx(i)) .* (y - my(i))); sxx = sxx + sum((x - mx(i)).^2);
    end
    beta = sxy / sxx; alpha = mean(my - beta * mx);
end

function rr = indiv_r(P, S, shift, minPts)
    ids = intersect(unique(P.id), unique(S.id)); rr = nan(numel(ids), 1);
    for i = 1:numel(ids)
        rr(i) = indiv_curve(P(P.id == ids(i), :), S(S.id == ids(i), :), shift, minPts);
    end
end

function rr = indiv_curve(Pi, Si, shifts, minPts)
    % Pearson correlation within one participant for each shift.
    rr = nan(size(shifts));
    [tp, ia] = unique(Pi.time);
    if numel(tp) < 2, return; end
    yp = log10(Pi.conc(ia));
    for j = 1:numel(shifts)
        tq = max(Si.time - shifts(j), tp(1)); ok = tq <= tp(end);
        if sum(ok) < minPts, continue; end
        a = log10(Si.conc(ok)); b = interp1(tp, yp, tq(ok));
        if std(a) == 0 || std(b) == 0, continue; end
        c = corrcoef(a, b); rr(j) = c(1, 2);
    end
end

function p = sign_test(nPos, nNeg)
    % Exact two-sided sign test (zeros excluded), no toolbox needed.
    n = nPos + nNeg; p = NaN;
    if n == 0, return; end
    k = 0:min(nPos, nNeg);
    p = min(1, 2 * sum(exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) - n*log(2))));
end

function Y = pct(X, p, dim)
    % Percentiles with linear interpolation, NaNs removed (vectors only for dim 1).
    X = X(:); X = sort(X(isfinite(X)));
    if isempty(X), Y = nan(1, numel(p)); return; end
    n = numel(X); Y = zeros(1, numel(p));
    for j = 1:numel(p)
        pos = 1 + (n - 1) * p(j) / 100; lo = floor(pos); hi = ceil(pos); w = pos - lo;
        Y(j) = X(lo) + w * (X(hi) - X(lo));
    end
    if nargin > 2 && dim ~= 1, Y = Y'; end
end

function v = quant(x, p)
    v = pct(x, 100*p, 1);
end

% ---- broken time axis ------------------------------------------------------------
function x = bx(th, B)
    % Map hours to the broken axis: [.. t1] unchanged, [t2 ..] moved left so
    % that it starts gap units after t1; values inside the gap become NaN
    % (lines break there).
    x = nan(size(th));
    m1 = th <= B.t1; x(m1) = th(m1);
    m2 = th >= B.t2; x(m2) = B.t1 + B.gap + (th(m2) - B.t2);
end

function bx_axis(ax, opt, xMin, ticks)
    if nargin < 4, ticks = opt.xTick; end
    set(ax, 'XLim', [xMin, bx(24.3, opt.B)], 'XTick', bx(ticks, opt.B), ...
        'XTickLabel', compose('%g', ticks), 'XTickLabelRotation', 0);
end

function bx_break(ax, B)
    % Break mark on the time axis (call after the y-limits are final).
    yl = ylim(ax); xm = B.t1 + B.gap/2;
    text(ax, xm, yl(1), '//', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'w', 'Margin', 0.1, 'FontSize', hs.FS_NOTE, 'Color', hs.C.ink2, 'Clipping', 'off');
end

function model_band(ax, F, k, which, mode, col, B)
    [~, lo, hi] = model_curve(F, k, which, mode);
    if isempty(lo), return; end
    th = F.unc.t / 60;
    for seg = {th <= B.t1, th >= B.t2}
        m = seg{1};
        if sum(m) >= 2, hs.band(ax, bx(th(m), B), lo(m), hi(m), col); end
    end
end

function model_line(ax, F, k, which, mode, col, ls, B)
    y = model_curve(F, k, which, mode);
    plot(ax, bx(F.sim.t/60, B), y, ls, 'Color', col, 'LineWidth', hs.LW_DATA);
end

function [y, lo, hi] = model_curve(F, k, which, mode)
    which = char(which); f = char(F.cytoNames(k));
    b = F.p.([which 'Baseline']).(f);
    y = F.sim.(which)(:, k); if mode == "fold", y = y / b; end
    lo = []; hi = [];
    if isempty(F.unc), return; end
    if mode == "abs"
        lo = F.unc.([which 'Lo'])(:, k); hi = F.unc.([which 'Hi'])(:, k);
    else
        lo = F.unc.([which 'FoldLo'])(:, k); hi = F.unc.([which 'FoldHi'])(:, k);
    end
end

function point_data(ax, S, col, mk, scale, B)
    % Means with CI (no line), drawn last.
    if isempty(S.t), return; end
    x = bx(S.t/60, B); ok = isfinite(x);
    if mk == 's', ec = col; ms = 3; else, ec = 'w'; ms = 4; end   % dense plasma squares: no white rim
    errorbar(ax, x(ok), S.m(ok)/scale, (S.m(ok) - S.lo(ok))/scale, (S.hi(ok) - S.m(ok))/scale, mk, ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', ec, 'MarkerSize', ms, 'CapSize', 0, ...
        'LineWidth', 0.8, 'LineStyle', 'none');
end

function line_data(ax, S, col, mk, dxMin, B)
    % Means with CI joined in time order (raw-data figure), broken at the gap.
    x = bx((S.t + dxMin)/60, B); ok = isfinite(x);
    m1 = S.t <= 480;                                   % join 0-8 h only, not across the gap
    plot(ax, x(m1), S.m(m1), '-', 'Color', col, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    errorbar(ax, x(ok), S.m(ok), S.m(ok) - S.lo(ok), S.hi(ok) - S.m(ok), mk, 'Color', col, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 2.5, 'CapSize', 0, ...
        'LineWidth', 0.5, 'LineStyle', 'none');
end

function draw_subjects(ax, S, col, B)
    for k = 1:numel(S.subj)
        d = S.subj{k};
        plot(ax, bx(d(:, 1)/60, B), d(:, 2), '-', 'Color', hs.tint(col, 0.80), 'LineWidth', 0.4, ...
            'HandleVisibility', 'off');
    end
end

function draw_lloq(ax, C, nm, ~)
    if isfield(C, 'lloqS') && isfield(C.lloqS, char(nm)) && isfinite(C.lloqS.(char(nm)))
        yline(ax, C.lloqS.(char(nm)), ':', 'Color', hs.C.ref, 'LineWidth', hs.LW_THIN, 'HandleVisibility', 'off');
    end
end

function finish_pair(axP, axS, k, opt, xMin)
    hs.axes(axP); hs.axes(axS);
    bx_axis(axP, opt, xMin, opt.xTickGrid); bx_axis(axS, opt, xMin, opt.xTickGrid);
    axP.XTickLabel = [];
    if mod(k - 1, opt.nCols) == 0
        ylabel(axP, 'Plasma (pg/mL)'); ylabel(axS, 'Saliva (pg/mL)');
    end
    bx_break(axP, opt.B); bx_break(axS, opt.B);
end

% ---- small helpers ----------------------------------------------------------------
function [iP, iS] = grid_index(c, nCols)
    block = ceil(c / nCols); col = mod(c - 1, nCols) + 1;
    iP = 2*(block - 1)*nCols + col; iS = iP + nCols;
end

function lim = log_lim(v, pad)
    v = v(isfinite(v) & v > 0);
    if isempty(v), lim = [1 10]; return; end
    lim = [min(v)/pad, max(v)*pad];
    if lim(2) <= lim(1), lim = lim .* [0.5 2]; end
end

function decade_ticks(ax)
    % Log axes: one tick per decade (every second decade if > 3), horizontal labels.
    for a = 'XY'
        if ~strcmp(ax.([a 'Scale']), 'log'), continue; end
        lim = ax.([a 'Lim']); e = ceil(log10(lim(1))):floor(log10(lim(2)));
        if numel(e) > 2, e = e(mod(e - e(end), 2) == 0); end
        if numel(e) >= 2
            ax.([a 'Tick']) = 10.^e; ax.([a 'TickLabel']) = compose('%g', 10.^e);
        end
    end
    ax.XTickLabelRotation = 0;
end

function lim = robust_lim(x)
    v = sort(log10(x(x > 0))); n = numel(v);
    lo = v(max(1, ceil(0.02*n))); hi = v(max(1, floor(0.98*n)));
    pad = 0.12*(hi - lo) + 0.1;
    lim = 10.^[lo - pad, hi + pad];
end

function c = arm_color(name)
    c = hs.C.ink2;
    for cand = [string(name), string(name) + "e"]
        try, c = hs.sem(char(cand)); return; catch, end
    end
end

function s = cyto_label(nm)
    m = containers.Map({'TNF','IL6','IL8','IL10','IL1RA','GCSF','IP10','MCP1','MIP1a','MIP1B'}, ...
        {'TNF','IL-6','IL-8','IL-10','IL-1RA','G-CSF','IP-10','MCP-1','MIP-1\alpha','MIP-1\beta'});
    key = char(nm); if isKey(m, key), s = m(key); else, s = key; end
end

function s = fmt_conc(x)
    if x >= 100, s = sprintf('%.0f', x); elseif x >= 10, s = sprintf('%.1f', x); else, s = sprintf('%.2g', x); end
end

function two_line_title(ax, main, sub)
    c = hs.C.ink2;
    title(ax, {main, sprintf('\\rm\\fontsize{%g}\\color[rgb]{%.3f,%.3f,%.3f}%s', hs.FS_NOTE, c(1), c(2), c(3), sub)});
end

function grid_legend(tl, items)
    ax = nexttile(tl, 1); hold(ax, 'on');
    h = gobjects(size(items, 1), 1);
    for i = 1:size(items, 1)
        h(i) = plot(ax, NaN, NaN, 'LineStyle', items{i, 4}, 'Marker', items{i, 3}, 'Color', items{i, 2}, ...
            'MarkerFaceColor', items{i, 2}, 'MarkerEdgeColor', 'w', 'LineWidth', hs.LW_DATA, 'MarkerSize', 4.5);
    end
    lgd = legend(ax, h, items(:, 1)', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', hs.FS_LABEL, 'AutoUpdate', 'off');
    lgd.Layout.Tile = 'north';
end

function h = corner_note(ax, str, loc, col)
    switch loc
        case 'nw', pos = [0.03 0.97]; ha = 'left';  va = 'top';
        case 'se', pos = [0.98 0.04]; ha = 'right'; va = 'bottom';
        otherwise, pos = [0.98 0.97]; ha = 'right'; va = 'top';
    end
    h = text(ax, pos(1), pos(2), str, 'Units', 'normalized', 'HorizontalAlignment', ha, ...
        'VerticalAlignment', va, 'FontSize', hs.FS_NOTE, 'FontName', hs.FONT, 'Color', col);
end
