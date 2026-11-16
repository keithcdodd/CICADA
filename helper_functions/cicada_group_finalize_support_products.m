function cicada_group_finalize_support_products(stage_root, support_root, move_fn)
% CICADA_GROUP_FINALIZE_SUPPORT_PRODUCTS Replace regenerated Group support tree.
% Successful reruns REPLACE prior products. A failed promotion attempts to
% restore the previous tree; any un-restorable backup is left on disk.
%
% move_fn is an optional function handle for isolated failure-injection
% tests. Production calls use MATLAB's normal movefile implementation.

if nargin < 3 || isempty(move_fn)
    move_fn = @movefile;
end
stage_root = char(string(stage_root));
support_root = char(string(support_root));
if strcmp(stage_root, support_root) || ...
        ~strcmp(fileparts(stage_root), fileparts(support_root))
    error('CICADA:SupportFinalizePaths', ...
        'Staging and canonical support folders must be distinct siblings.');
end
if ~isfolder(stage_root)
    error('CICADA:SupportStageMissing', ...
        'Staged support products are missing: %s', stage_root);
end
if isfile(support_root)
    error('CICADA:SupportTargetIsFile', ...
        'Cannot replace a regular file with support products: %s', support_root);
end

backup_root = '';
if isfolder(support_root)
    backup_root = [tempname(fileparts(support_root)), '_previous_support'];
    [ok, detail] = local_try_move(move_fn, support_root, backup_root);
    if ~ok
        error('CICADA:SupportBackupMoveFailure', ...
            ['Cannot preserve the previous support products: %s. ' ...
             'New staged products were not installed.'], detail);
    end
end

[ok, detail] = local_try_move(move_fn, stage_root, support_root);
if ~ok
    % A failing move may leave a partially created canonical destination.
    % Quarantine it rather than deleting anything needed for investigation.
    partial_root = '';
    if isfolder(support_root)
        partial_root = [tempname(fileparts(support_root)), '_incomplete_support'];
        [partial_ok, partial_detail] = local_try_move( ...
            move_fn, support_root, partial_root);
        if ~partial_ok
            error('CICADA:SupportRestoreBlocked', ...
                ['New support promotion failed: %s. Incomplete target at %s ' ...
                 'could not be moved aside: %s. Previous output backup: %s. ' ...
                 'Manual recovery is required.'], ...
                detail, support_root, partial_detail, backup_root);
        end
    end
    if ~isempty(backup_root)
        [restore_ok, restore_detail] = local_try_move( ...
            move_fn, backup_root, support_root);
        if ~restore_ok
            error('CICADA:SupportRestoreFailure', ...
                ['New support promotion failed: %s. Restoration failed: %s. ' ...
                 'Previous output remains at %s. Incomplete output at %s.'], ...
                detail, restore_detail, backup_root, partial_root);
        end
    end
    error('CICADA:SupportOutputMoveFailure', ...
        ['Could not finalize Group support products: %s. ' ...
         'Prior outputs were restored when present. Check staging ' ...
         'at %s. Any incomplete output was preserved at %s.'], ...
        detail, stage_root, partial_root);
end

% Only the old, successfully replaced canonical tree is deleted. If that
% cleanup fails, keep the new result and warn about the leftover backup.
if ~isempty(backup_root) && isfolder(backup_root)
    try
        [deleted, reason] = rmdir(backup_root, 's');
    catch ME
        deleted = false;
        reason = ME.message;
    end
    if ~deleted
        warning('CICADA:SupportBackupCleanupFailure', ...
            'New support outputs finalized, but old backup remains at %s: %s', ...
            backup_root, char(string(reason)));
    end
end
end

function [ok, detail] = local_try_move(move_fn, source, destination)
try
    [ok, detail] = move_fn(source, destination);
    ok = logical(ok);
    detail = char(string(detail));
catch ME
    ok = false;
    detail = ME.message;
end
end
