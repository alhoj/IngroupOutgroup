function calculate_twISC(cfg)
% Calculates individual ISC maps with respect to a reference group
% 
% Usage:
%   calculate_twISC(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects of interest
% 
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
%
%   cfg.toi = time points of interest in a vector
% 
%   cft.outdir = label of the output folder where the individual ISC nifitis are saved


%% Input validation

if ~ismember(cfg.subs,{'heteros','homos','both'})
    error('cfg.subs should be either ''heteros'', ''homos'' or ''both''!')
end
if ~ismember(cfg.nii,{'epi_movie.nii','epi_movie1.nii','epi_movie2.nii'})
    error('cfg.nii should be either ''epi_movie.nii'', ''epi_movie1.nii'' or ''epi_movie2.nii''!')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['Could not find mask: ' cfg.mask])
end
if ~isvector(cfg.toi)
    error('cfg.toi must be a vector of numbers, e.g. [1 2 3 4 5]!')
end

cfg.toi=sort(cfg.toi);

disp(['Calculating time windowed ISC for ' cfg.subs])
disp(['Using time points: ' num2str(cfg.toi)])
fprintf('\n')
%% Read IDs for the subject groups

codes=importdata('/m/nbe/scratch/braindata/jaalho/gaypriest/subIDs.txt'); % Import the subjects code text to split them into patients, controls and reference
% Split codes into patient codes, control codes and reference codes

mode=0; %  1 for heteros, 2 for homos
sample=0; % Increasing index for each subject in the same group
for codei=1:length(codes)
    
    if ~strcmp(codes{codei}(1:3),'sub') % If it is not a subject code
        codes{codei};
        sample=0; % Set the indexing number to zero
        mode=mode+1; % Increase the mode by 1 to get to the next category
    else
        sample=sample+1;  % Increase the index of the subject in this category
        if mode==1
            heteros{sample}=codes{codei};
        elseif mode==2
            homos{sample}=codes{codei};
        end
    end
end

% Set IDs for subjects of interest and reference subjects
switch cfg.subs
    case 'heteros'
        cfg.infiles=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii);
    case 'homos'
        cfg.infiles=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii);
    case 'both'
        cfg.infiles=[strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii) strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii)];
end


%% Load mask
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/');

disp('Loading mask...')
mask=load_nii(cfg.mask);
inmask=find(mask.img);
fprintf('\n')

%% Load brain data

% Load an example NIFTI file to find the dimensions needed for
% preallocating the matrices (this will speed things up)
test_nii=load_nii(cfg.infiles{1});

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects
ntps=size(temp(:,inmask),1);

assert(cfg.toi(end)<=ntps)
nvox=size(temp(:,inmask),2);
nsub=length(cfg.infiles);
ntoi=length(cfg.toi);

% Preallocate data matrices with the dimensions of time points of interest x subjects x voxels 
allsubs=zeros(ntoi,nsub,nvox); 

disp('Loading brain data...')
for i=1:nsub
    disp([num2str(i) ' out of ' num2str(nsub)])
    nii=load_nii(cfg.infiles{i});
    temp=permute(nii.img,[4 1 2 3]);
    allsubs(:,i,:)=zscore(temp(cfg.toi,inmask));
%     allsubs(:,i,:)=temp(:,inmask);
end
fprintf('\n')

%% Calculate ISC

% Preallocate the correlation matrix
cormat=zeros(nvox,nsub,nsub);

disp('Calculating ISC...')
for voxi=1:nvox
    if mod(voxi,1000)==0 % Show the status every 1000 voxels
        disp([num2str(voxi) '/' num2str(nvox) ' voxels'])
    end
       
    % Calculate correlation over time
    cormat(voxi,:,:)=corr(squeeze(allsubs(:,:,voxi)));
  
end
fprintf('\n')

%% Average to obtain indivual ISC maps and save them
disp('Saving correlation matrix...')

% Check if output directory exists; if not, create it
if ~exist(cfg.outdir,'dir') 
    system(['mkdir -p ' cfg.outdir]);
end

% Replace possible NaN values with zeros.
cormat(isnan(cormat))=0;

filename=[cfg.outdir '/cormat.mat'];
save(filename,'cormat');

disp('Done!')
fprintf('\n')