function save_figure_robust(fig, filename, dpi)
%SAVE_FIGURE_ROBUST Save figure reliably across platforms and MATLAB versions.
%
%   save_figure_robust(fig, filename, dpi)
%
%   Inputs:
%     - fig:      Handle to a figure (required)
%     - filename: Output file path (e.g., 'myplot.png', 'fig.pdf') (required)
%     - dpi:      Output resolution in dots per inch (default = 300)
%
%   This function uses 'print' instead of 'saveas', configures appropriate
%   renderer, and works in headless/server environments.

    if nargin < 3 || isempty(dpi)
        dpi = 300;
    end

    % Ensure software OpenGL to avoid GPU/display issues on servers
    try
        opengl('save', 'software');
        opengl('software');
    catch
        warning('Could not set OpenGL to software mode.');
    end

    % Enforce consistent appearance
    set(fig, 'Visible', 'off', ...
             'PaperPositionMode', 'auto', ...
             'InvertHardcopy', 'off', ...
             'Color', 'w');  % White background

    % Standardize font size (adjust as needed)
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12);

    % Determine output format
    [~, ~, ext] = fileparts(filename);
    ext = lower(ext);

    % Set renderer based on output type
    switch ext
        case {'.png', '.jpg', '.jpeg', '.tiff', '.bmp'}
            set(fig, 'Renderer', 'opengl');
            % Determine output format
			[~, ~, ext] = fileparts(filename);
			ext = lower(ext);

			switch ext
    			case '.png'
        			set(fig, 'Renderer', 'opengl');
        			fmt = '-dpng';
    			case {'.jpg', '.jpeg'}
        			set(fig, 'Renderer', 'opengl');
        			fmt = '-djpeg';
    			case '.tiff'
        			set(fig, 'Renderer', 'opengl');
        			fmt = '-dtiff';
    			case '.bmp'
        			set(fig, 'Renderer', 'opengl');
        			fmt = '-dbmp';
    			case '.pdf'
        			set(fig, 'Renderer', 'painters');
        			fmt = '-dpdf';
    			case '.eps'
        			set(fig, 'Renderer', 'painters');
        			fmt = '-depsc';
    			otherwise
        			warning('Unknown extension "%s". Defaulting to PNG.', ext);
        			set(fig, 'Renderer', 'opengl');
        			fmt = '-dpng';
        			filename = [filename '.png'];
			end
        case {'.pdf'}
            set(fig, 'Renderer', 'painters');  % Use vector graphics
            fmt = '-dpdf';
        case {'.eps'}
            set(fig, 'Renderer', 'painters');  % Use vector graphics
            fmt = '-depsc';
        otherwise
            warning('Unknown extension "%s". Defaulting to PNG.', ext);
            set(fig, 'Renderer', 'opengl');
            fmt = '-dpng';
            filename = [filename '.png'];
    end

    % Save figure using print
    try
        print(fig, filename, fmt, ['-r' num2str(dpi)]);
    catch ME
        error('Failed to save figure: %s', ME.message);
    end

    % Optional: close figure after saving
    % close(fig);
end
