function compare_file = find_compare_file(output_dir, compare_file, valid_tags)
% FIND_COMPARE_FILE Resolve a comparison file or comparison tag.
%
% compare_file may be:
%   1. A full path to an existing comparison file, or
%   2. A valid comparison tag such as '8p', '12p', etc.
%
% If a valid tag is supplied and the corresponding cleaned file already
% exists, return its full path.
%
% If a valid tag is supplied but the file does NOT yet exist, preserve and
% return the tag so that a caller such as Manual_CICADA can recreate it.
%
% Invalid/missing specifications fall back to the standard '8p' comparison.

cleaned_dir = fullfile(output_dir, 'cleaned');

if nargin < 3 || isempty(valid_tags)
    valid_tags = {'6p', '8p', '9p', '12p', '16p', '18p', ...
        '24p', '28p', '30p', '32p', '36p'};
end

if nargin < 2 || isempty(compare_file) || ~ischar(compare_file)

    fprintf('No valid compare file/tag supplied. Using standard 8p.\n')
    compare_file = '8p';

end


% If a full existing filepath was supplied, use it directly.
if isfile(compare_file)
    return
end


% Otherwise it must be a supported comparison tag.
if ~ismember(compare_file, valid_tags)

    fprintf(['Compare specification is not an existing file or valid tag. ' ...
        'Using standard 8p.\n'])

    compare_file = '8p';

end

compare_tag = compare_file;


% See whether this tagged comparison already exists.
compare_file_info = ...
    dir(fullfile(cleaned_dir, ['*_', compare_tag, '_*.nii.gz']));

if isscalar(compare_file_info)

    compare_file = fullfile( ...
        compare_file_info(1).folder, ...
        compare_file_info(1).name);

elseif numel(compare_file_info) > 1

    error('CICADA:AmbiguousCompareFile', ...
        ['Multiple cleaned comparison files were found for tag "%s" in:\n%s\n' ...
         'Please provide the intended full filepath explicitly.'], ...
        compare_tag, cleaned_dir);

else

    % IMPORTANT: preserve the valid tag. A downstream caller may be able
    % to recreate the comparison from the saved regression design matrix.
    fprintf(['Comparison file for %s does not currently exist. ' ...
        'Returning the valid tag so it can be recreated if needed.\n'], ...
        compare_tag);

    compare_file = compare_tag;

end

end