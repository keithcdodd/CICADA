function files = cicada_find_files(pattern)
%CICADA_FIND_FILES List visible regular files matching a DIR pattern.
%   FILES = CICADA_FIND_FILES(PATTERN) uses MATLAB's DIR wildcard matching,
%   then excludes directories and files whose names begin with a dot.
%   This avoids accidental matches to macOS AppleDouble (._*) files.
%
%   PATTERN can be a relative or absolute path with DIR wildcards, e.g.:
%       cicada_find_files('*Auto.*')
%       cicada_find_files(fullfile(cleaned_dir, '*CICADA*auto*.nii.gz'))
%
%   Returns a DIR-compatible struct array (possibly empty). This function
%   does not delete or move anything, and it does not search recursively
%   unless PATTERN explicitly uses DIR's recursive ** syntax.
%
    files = dir(pattern);
    if isempty(files)
        return
    end

    is_visible_file = ~startsWith({files.name}, '.') & ~[files.isdir];
    files = files(is_visible_file);
end
