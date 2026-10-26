function info = cicada_git_provenance(source_file)
%CICADA_GIT_PROVENANCE Identify the Git revision of the code being run.
%   INFO = CICADA_GIT_PROVENANCE(SOURCE_FILE) examines the Git checkout
%   containing SOURCE_FILE (e.g. mfilename('fullpath') in AutoCICADA).
%   If SOURCE_FILE is omitted, the helper's own location is used.
%   No repository files are changed. Missing Git metadata never errors.
%
    info = struct('Available', false, ...
                  'Status', 'unavailable');

    if nargin < 1 || isempty(source_file)
        source_file = mfilename('fullpath');
    end

    source_dir = fileparts(char(source_file));
    if ~isfolder(source_dir)
        info.Status = 'source_directory_not_found';
        return
    end

    previous_dir = pwd;
    restore_dir = onCleanup(@() cd(previous_dir)); %#ok<NASGU>
    cd(source_dir);

    [status, root] = system('git rev-parse --show-toplevel 2>/dev/null');
    if status ~= 0 || isempty(strtrim(root))
        info.Status = 'not_a_git_checkout_or_git_unavailable';
        return
    end

    [status, commit] = system('git rev-parse --verify HEAD 2>/dev/null');
    if status ~= 0 || isempty(strtrim(commit))
        info.Status = 'git_commit_unavailable';
        return
    end

    info.Available = true;
    info.Status = 'ok';
    info.Commit = strtrim(commit);

    [status, branch] = system( ...
        'git symbolic-ref --quiet --short HEAD 2>/dev/null');
    if status == 0 && ~isempty(strtrim(branch))
        info.Branch = strtrim(branch);
    else
        info.Branch = 'detached_HEAD';
    end

    [status, tag] = system( ...
        'git describe --tags --exact-match HEAD 2>/dev/null');
    if status == 0 && ~isempty(strtrim(tag))
        info.Tag = strtrim(tag);
    end

    % Include tracked modifications and ordinary untracked files.
    % A dirty checkout cannot be reconstructed from its commit alone.
    [status, changes] = system( ...
        'git status --porcelain --untracked-files=normal 2>/dev/null');
    if status == 0
        info.Dirty = ~isempty(strtrim(changes));
    else
        info.Status = 'git_status_unavailable';
    end
end
