function save_figure_robust(fig, filename, dpi)
%SAVE_FIGURE_ROBUST Save a figure reliably across platforms and MATLAB versions.
%
%   save_figure_robust(fig, filename, dpi)
%
%   Inputs:
%     - fig:      Handle to a figure (required)
%     - filename: Output file path (e.g., 'myplot.png', 'fig.pdf') (required)
%     - dpi:      Output resolution in dots per inch (default = 300)
%
%   Features:
%     - Automatically selects correct renderer (OpenGL or Painters)
%     - Handles missing extensions by defaulting to PNG
%     - Works on headless servers via software OpenGL
%     - Applies consistent formatting (font size, white background)

    if nargin < 3 || isempty(dpi)
        dpi = 300;
    end

    % Ensure figure handle is valid
    if ~ishandle(fig) || ~strcmp(get(fig, 'Type'), 'figure')
        error('First argument must be a valid figure handle.');
    end

   try
        if usejava('desktop') && feature('opengl','supported')
            opengl('software');
        else
            warning('OpenGL not supported or no desktop environment; skipping OpenGL setup. This is usually not an issue. Check QC plots!');
        end
    catch
        warning('Could not set OpenGL to software mode. This is usually not an issue. Check QC plots.');
    end

    % Standard appearance
    set(fig, 'Visible', 'off', ...
             'PaperPositionMode', 'auto', ...
             'InvertHardcopy', 'off', ...
             'Color', 'w');  % white background
    %set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);

    % Handle missing extension (default to PNG)
    [filepath, name, ext] = fileparts(filename);
    if isempty(ext)
        ext = '.png';
        filename = fullfile(filepath, [name ext]);
    end
    ext = lower(ext);

    % Select format and renderer
    switch ext
        case '.png'
            set(fig, 'Renderer', 'opengl');  fmt = '-dpng';
        case {'.jpg', '.jpeg'}
            set(fig, 'Renderer', 'opengl');  fmt = '-djpeg';
        case '.tiff'
            set(fig, 'Renderer', 'opengl');  fmt = '-dtiff';
        case '.bmp'
            set(fig, 'Renderer', 'opengl');  fmt = '-dbmp';
        case '.pdf'
            set(fig, 'Renderer', 'painters'); fmt = '-dpdf';
        case '.eps'
            set(fig, 'Renderer', 'painters'); fmt = '-depsc';
        case '.svg'
            set(fig, 'Renderer', 'painters'); fmt = '-dsvg';  % MATLAB R2020a+
        otherwise
            warning('Unknown extension "%s". Defaulting to PNG.', ext);
            set(fig, 'Renderer', 'opengl');
            fmt = '-dpng';
            filename = fullfile(filepath, [name '.png']);
    end

    % Save the figure
    try
        print(fig, filename, fmt, ['-r' num2str(dpi)]);
        fprintf('Figure saved to: %s\n', filename);
    catch ME
        error('Failed to save figure: %s', ME.message);
    end

    % Optional: close figure automatically after saving
    % close(fig);
end
