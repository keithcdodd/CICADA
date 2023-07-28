% CADICA matlab only Denoising idea

% potential theory: (1) Run an ICA (e.g., reconstructive ICA) on all the voxel
% time series only falling within the gray matter mask. (2) We then
% correlate each IC to all the voxel timeseries. (3) We mark noise as any
% correlations where the top 5% are outside of the gray matter mask (or
% some variation on that theme... e.g., we could have a definitely outside
% of gray matter mask, and if any of that top 5% are within that, we remove
% it). OR (3.A) For each IC, we measure Grey Matter vs not Grey Matter
% overlap, and then also bandpass 0.008 to 0.15 Hz and compare. If Grey
% Matter proportion decreases, that IC is noise.
% Regardless, GM proportion must be larger than out of GM proportion. We
% could also take each IC and bandpass the ICs AND the functional file
% before performing correlation to check overlap (would help focus on the
% true source). Either way, we would expect GM overlap to significantly INCREASE if it is
% a true BOLD signal source in that case. (4) We run a nonaggressive
% regression from there (there should be a way to do this in matlab..., but
% if not, we can use FSL). (5) Then bandpass filtering and smoothing
% (ideally we have a way to apply these, with the regression, all at once, 
% or at least the bandpass filtering and smoothing at the same time)


% let's say we want to instead just load the nifti and run analysis from
% there: 
% gunzip(file.nii.gz)
% file = niftiread('filename.nii')
% gzip(file.nii.gz)
% inbrainmask_1D = sum(file>0,4) == 300;
% inbrainmask_4D = repmat(inbrainmask_1D, 1, 1, 1, 300);
% GMmask_4D, GMmask_1D, OutOfGMmask_4D, OutofGMmask_1D

clearvars

%%%%%%%%% set up that user may need to be adjust %%%%%%%%%%%%%%%%%%%
cadicafol = '/home/keithdodd/ExampleDataLocal/CADICA_Updated';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '108' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'102'};
sessions = {'01'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ultimately, want ICs to lean towards GM overlap instead of notGM, and we
% want this overlap to improve after bandpass filtering. The peaks/sources
% should like in GM, or at least not in NOT GM

addpath(cadicafol)

for j = 1:length(subjects)
    cd(cadicafol)
    currsubjfol = ['sub-', subjects{j}];
    cd(currsubjfol)
    for k = 1:length(sessions)
       cd(cadicafol)
       cd(currsubjfol)
       currsessfol = ['ses-', sessions{k}];
       cd(currsessfol)

       % grab functional file and functional mask
       

       % grab GM probability file (would need to be resampled)
       % This is easiest by calling afni or fsl, e.g., 3dresample. This
       % should definitely be possible here (e.g., calling a bash command)
       

       % Once you have resampled GM probability, you can mask at >0.67 to
       % get GM mask, and also mask at < 0.33 to get an OutofGM mask
       % (multiply both by functional mask so we stay within there)
       

       % Now, grab functional file within GMMaskm and run ICA
       % (reconstructive ICA should work, but there also used to be a
       % FastICA function add on, which might be very nice simply for 
       % estimating number of components). Could call MELODIC in fsl, and
       % then tell it to run in GM mask only... just need to see how easy
       % it is to grab timeseries from there for your own regression.
       

       % Run regression of ICs to the full functional file, as well as
       % regression of bandpassed ICs to bandpassed functional file


       % Mask for the top 5% of correlation values, mask for GM and
       % outofGM and sum for both. Do this for the regular, and bandpassed
       % versions. 


       % Check if (1) Bandpassing helped get more GM overlap, and (2) if
       % the bandpassed version has more GM overlap than not GM overlap.
       % All others, mark as noise.


       % Find/run some sort of nonaggressive regression. May be easiest to
       % call fsl for fsl_regfilt.





    end
end