function save_figure_robust(fig, filename, dpi)
%SAVE_FIGURE_ROBUST Save figure robustly for MATLAB R2020a+.
%
% Strategy:
% - exportgraphics() first
% - fallback to print()
% - no dependence on opengl/renderer switching
% - keeps figure hidden

    if nargin < 3 || isempty(dpi)
        dpi = 300;
    end

    if ~ishandle(fig) || ~strcmp(get(fig, 'Type'), 'figure')
        error('First argument must be a valid figure handle.');
    end

    set(fig, 'Visible', 'off', ...
             'Color', 'w', ...
             'InvertHardcopy', 'off', ...
             'PaperPositionMode', 'auto');

    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        try
            set(ax(k), 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        catch
        end
    end

    drawnow;

    [filepath, name, ext] = fileparts(filename);
    if isempty(ext)
        ext = '.png';
        filename = fullfile(filepath, [name ext]);
    end
    ext = lower(ext);

    try
        switch ext
            case '.png'
                try
                    exportgraphics(fig, filename, 'Resolution', dpi);
                catch
                    print(fig, filename, '-dpng', ['-r' num2str(dpi)]);
                end

            case {'.jpg', '.jpeg'}
                try
                    exportgraphics(fig, filename, 'Resolution', dpi);
                catch
                    print(fig, filename, '-djpeg', ['-r' num2str(dpi)]);
                end

            case {'.tif', '.tiff'}
                try
                    exportgraphics(fig, filename, 'Resolution', dpi);
                catch
                    print(fig, filename, '-dtiff', ['-r' num2str(dpi)]);
                end

            case '.pdf'
                try
                    exportgraphics(fig, filename, 'ContentType', 'vector');
                catch
                    print(fig, filename, '-dpdf');
                end

            case '.eps'
                print(fig, filename, '-depsc');

            case '.svg'
                try
                    exportgraphics(fig, filename, 'ContentType', 'vector');
                catch ME_svg
                    error('SVG export failed for "%s": %s', filename, ME_svg.message);
                end

            case '.bmp'
                print(fig, filename, '-dbmp', ['-r' num2str(dpi)]);

            otherwise
                warning('Unknown extension "%s". Defaulting to PNG.', ext);
                filename = fullfile(filepath, [name '.png']);
                try
                    exportgraphics(fig, filename, 'Resolution', dpi);
                catch
                    print(fig, filename, '-dpng', ['-r' num2str(dpi)]);
                end
        end

    catch ME
        error('Failed to save figure "%s": %s', filename, ME.message);
    end
end