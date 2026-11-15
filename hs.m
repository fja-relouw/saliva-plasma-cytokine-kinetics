classdef hs
% HS  Figure house style and export in the journal format.
%
%   fig = hs.figure(174, 0.56);      % width in mm, aspect = height/width
%   hs.axes(ax);                     % axis style (after plotting)
%   hs.panelLabels([ax1 ax2]);       % bold A, B
%   hs.save(fig, 'results/figures/Figure2');
%
% hs.save writes vector PDF, SVG and EPS and 600-dpi PNG and LZW-compressed TIFF at
% print size, in the journal style: Arial, no boxes, no titles above the artwork,
% midline decimal point (0·5) and true minus sign, transparency flattened (for EPS).

    properties (Constant)
        FONT     = 'Source Sans 3'      % drawing font; exported figures use Arial
        FS_PANEL = 9
        FS_LABEL = 8
        FS_TICK  = 7
        FS_NOTE  = 7
        LW_DATA  = 1.25
        LW_THIN  = 0.6
        LW_AXIS  = 0.6
        MS       = 4
        C = struct('ink', [0.110 0.102 0.141], 'ink2', [0.353 0.333 0.388], ...
            'axis', [0.243 0.227 0.278], 'ref', [0.553 0.533 0.588], 'context', [0.733 0.714 0.757], ...
            'accent', [0.761 0.251 0.478], 'amber', [0.878 0.569 0.247], ...
            'violet', [0.357 0.227 0.576], 'twilight', [0.373 0.651 0.878])
    end

    methods (Static)
        function c = sem(name)
            % Colour per treatment group.
            switch lower(name)
                case 'placebo',       c = hs.C.ink;
                case 'anakinra',      c = hs.C.amber;
                case 'dexamethasone', c = hs.C.violet;
                case 'tocilizumab',   c = hs.C.twilight;
                otherwise, error('hs:sem', 'No colour for "%s".', name);
            end
        end

        function m = ordered(n)
            % n colours for an ordered variable (time): light lilac -> night.
            h = {'#D9D3EA', '#9E92C9', '#5B4F99', '#26214A'};
            P = cell2mat(cellfun(@(x) sscanf(x(2:end), '%2x%2x%2x', [1 3]) / 255, h(:), 'UniformOutput', false));
            m = interp1(linspace(0, 1, 4), P, linspace(0, 1, 256));
            m = m(round(linspace(0.12, 1, n) * 255) + 1, :);
        end

        function c = tint(c, amount)
            % Mix with white.
            c = c + (1 - c) * amount;
        end

        function fig = figure(w, aspect)
            sz = [w/10 w*aspect/10];                 % cm
            fig = figure('Units', 'centimeters', 'Color', 'w', 'Position', [2 2 sz], ...
                'PaperUnits', 'centimeters', 'PaperSize', sz, 'PaperPosition', [0 0 sz]);
            set(fig, 'DefaultAxesFontName', hs.FONT, 'DefaultTextFontName', hs.FONT, ...
                'DefaultAxesColorOrder', [hs.C.ink; hs.C.amber; hs.C.violet; hs.C.twilight], ...
                'DefaultLineLineWidth', hs.LW_DATA, 'DefaultLineMarkerSize', hs.MS, ...
                'DefaultAxesFontSize', hs.FS_TICK, 'DefaultTextFontSize', hs.FS_NOTE);
            set(fig, 'WindowStyle', 'normal', 'UserData', struct('hsSizeCm', sz));
        end

        function axes(ax)
            set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.015], ...
                'LineWidth', hs.LW_AXIS, 'XColor', hs.C.axis, 'YColor', hs.C.axis, ...
                'FontName', hs.FONT, 'FontSize', hs.FS_TICK, ...
                'LabelFontSizeMultiplier', hs.FS_LABEL / hs.FS_TICK, ...
                'TitleFontSizeMultiplier', hs.FS_LABEL / hs.FS_TICK, ...
                'TitleFontWeight', 'bold', 'TitleHorizontalAlignment', 'left', ...
                'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
                'Layer', 'top', 'Color', 'w');
            ax.XLabel.Color = hs.C.ink; ax.YLabel.Color = hs.C.ink;
            for r = [ax.XAxis(:); ax.YAxis(:)]'
                if strcmp(r.Scale, 'log'), hs.logTicks(r); end
            end
            if ~isempty(ax.Legend), set(ax.Legend, 'Box', 'off', 'FontSize', hs.FS_NOTE); end
        end

        function logTicks(r)
            % Plain-number ticks on a log axis: decades, or 1-2-5 / 1-3 steps for short ranges.
            lim = r.Limits;
            if any(~isfinite(lim)) || lim(1) <= 0, return; end
            d0 = floor(log10(lim(1))); d1 = ceil(log10(lim(2)));
            if d1 - d0 <= 1, base = [1 2 5]; elseif d1 - d0 <= 2, base = [1 3]; else, base = 1; end
            t = kron(10.^(d0:d1), base);
            t = t(t >= lim(1)*(1 - 1e-9) & t <= lim(2)*(1 + 1e-9));
            if numel(t) < 2, return; end
            lab = cell(size(t));
            for i = 1:numel(t)
                if t(i) >= 1e5 || t(i) < 1e-3, lab{i} = sprintf('10^{%d}', round(log10(t(i))));
                else, lab{i} = sprintf('%g', t(i)); end
            end
            r.TickValues = t; r.TickLabels = lab; r.MinorTick = 'on';
        end

        function h = band(ax, x, lo, hi, c)
            % Confidence band: opaque tint, no edge.
            h = fill(ax, [x(:); flipud(x(:))], [lo(:); flipud(hi(:))], hs.tint(c, 0.82), ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        function panelLabels(axs)
            L = 'ABCDEFGHIJKLMNOP';
            for k = 1:numel(axs)
                ax = axs(k); u = ax.Units; ax.Units = 'centimeters'; p = ax.Position; ax.Units = u;
                text(ax, -0.9 / p(3), 1, L(k), 'Units', 'normalized', 'FontSize', hs.FS_PANEL, ...
                    'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
                    'Color', hs.C.ink, 'Clipping', 'off');
            end
        end

        function save(fig, name)
            % Export at print size in the journal style (see the header).
            folder = fileparts(name);
            if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
            sz = fig.UserData.hsSizeCm;
            fig.WindowStyle = 'normal'; fig.Units = 'centimeters';
            p = fig.Position; fig.Position = [p(1:2) sz];
            set(fig, 'PaperUnits', 'centimeters', 'PaperSize', sz, 'PaperPosition', [0 0 sz]);
            drawnow; hs.journalStyle(fig); drawnow;
            exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
            exportgraphics(fig, [name '.eps'], 'ContentType', 'vector', 'BackgroundColor', 'white');
            try
                print(fig, [name '.svg'], '-dsvg', '-vector');      % R2022b and later
            catch
                print(fig, [name '.svg'], '-dsvg', '-painters');
            end
            for ext = {'.png', '.tif'}
                file = [name ext{1}];
                exportgraphics(fig, file, 'Resolution', 600, 'BackgroundColor', 'white');
                img = imread(file);                    % write the resolution into the file
                if strcmp(ext{1}, '.png')
                    imwrite(img, file, 'ResolutionUnit', 'meter', 'XResolution', round(600/0.0254), ...
                        'YResolution', round(600/0.0254));
                else
                    imwrite(img, file, 'Resolution', 600, 'Compression', 'lzw');
                end
            end
            fprintf('Saved %s\n', name);
        end

        function journalStyle(fig)
            set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
            for tl = findall(fig, 'Type', 'tiledlayout')'
                tl.Title.String = '';
            end
            delete(findall(fig, 'Type', 'subplottext'));
            for ax = findall(fig, 'Type', 'axes')'
                ax.Box = 'off';
                for r = [ax.XAxis(:); ax.YAxis(:)]'
                    if isempty(r.TickLabels), continue; end
                    v = r.TickValues; lab = r.TickLabels;
                    r.TickValues = v; r.TickLabels = hs.journalNumbers(lab);
                end
            end
            for lg = findall(fig, 'Type', 'legend')'
                lg.Box = 'off'; lg.String = hs.journalNumbers(lg.String);
            end
            for t = findall(fig, 'Type', 'text')'
                if ~strcmpi(t.Interpreter, 'latex'), t.String = hs.journalNumbers(t.String); end
            end
            for h = findall(fig, '-property', 'FaceAlpha')'   % flatten transparency, keep behind
                a = h.FaceAlpha;
                if isnumeric(a) && isscalar(a) && a < 1 && isnumeric(h.FaceColor) && numel(h.FaceColor) == 3
                    h.FaceColor = a * h.FaceColor + (1 - a) * [1 1 1]; h.FaceAlpha = 1;
                    try, uistack(h, 'bottom'); catch, end
                end
            end
        end

        function s = journalNumbers(s)
            % 0.5 -> 0·5 and -2 -> −2, leaving TeX commands with numeric arguments alone.
            if iscell(s), s = cellfun(@hs.journalNumbers, s, 'UniformOutput', false); return; end
            if isstring(s), s = string(hs.journalNumbers(cellstr(s))); return; end
            if ~ischar(s) || isempty(s), return; end
            if size(s, 1) > 1, s = char(hs.journalNumbers(cellstr(s))); return; end
            mid = char(183); minus = char(8722);
            [cmd, parts] = regexp(s, '\\(?:color\[rgb\]|fontsize|color)\{[^}]*\}', 'match', 'split');
            for i = 1:numel(parts)
                p = regexprep(parts{i}, '(?<=\d)\.(?=\d)', mid);
                parts{i} = regexprep(p, ['(^|[\s(\[{,=<>~])-(?=[\d\.' mid '])'], ['$1' minus]);
            end
            s = parts{1};
            for i = 1:numel(cmd), s = [s cmd{i} parts{i + 1}]; end %#ok<AGROW>
        end
    end
end
