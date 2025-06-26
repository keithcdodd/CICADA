function compare_file = find_compare_file(output_dir, compare_file, valid_tags)
% will look for compare_file, will handle if it is instead a valid
% compare_tag and find the corresponding file if it exists already! Otherwise just return it as a valid tag 
% valid tags like {'8p', '9p'} etc.
% point is to look for a compare file given the compare_tag, and if not
% there, try to make it. 
% compare_file needs to be a char array
% if it already is a full filepath that exists, that's great!
% only use for manual or group

if ~ischar(compare_file)
    fprintf('Compare file is not a char array. Making it the default...\n')
    compare_file = '';
end

cleaned_dir = [output_dir, '/cleaned'];
task_dir = output_dir;

if ~ismember(compare_file, valid_tags)
    if ~isfile(compare_file)
        fprintf('Will compare to standard 8 parameter and Auto CICADA \n')
        compare_file = ''; 
    end
else
    compare_tag = compare_file;
    fprintf(['Will try to compare to standard ', compare_tag, ', if it exists! \n'])

    % Auto CICADA should have made the necessary files to compute if
    % the file does not exist:
    compare_file_info = dir([cleaned_dir, '/*', compare_file, '*']); % help see if file exists as it should

    if isempty(compare_file_info)
        % that means the file does not currently exist
        compare_file = '';
    else
        compare_file = [compare_file_info.folder, '/', compare_file_info.name]; % Update it to match to the tagged type of comparison denoiser
    end
end

end