function calculate_BOLD(cfg)
% Calculates BOLD signals for heteros and homos


%% Input validation

if ~ismember(cfg.nii,{'epi_movie.nii','epi_movie1.nii','epi_movie2.nii'})
    error('cfg.nii should be either ''epi_movie.nii'', ''epi_movie1.nii'' or ''epi_movie2.nii''!')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['Could not find mask: ' cfg.mask])
end
if ~isvector(cfg.toi) && ~isempty(cfg.toi)
    error('cfg.toi must be a vector of numbers, e.g. [1 2 3 4 5]!')
end


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

cfg.infilesHt=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii);
cfg.infilesHo=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii);
    

%% Load mask
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/');

disp('Loading mask...')
mask=load_nii(cfg.mask);
inmask=find(mask.img);
fprintf('\n')

%% Load brain data

% Load an example NIFTI file to find the dimensions needed for
% preallocating the matrices (this will speed things up)
test_nii=load_nii(cfg.infilesHt{1});

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects
ntps=size(temp(:,inmask),1);

if isempty(cfg.toi)
    cfg.toi=1:ntps;
end

assert(cfg.toi(end)<=ntps)
nvox=size(temp(:,inmask),2);
nHt=length(cfg.infilesHt);
nHo=length(cfg.infilesHo);
ntoi=length(cfg.toi);

% Preallocate data matrices with the dimensions of time points of interest x subjects x voxels 
allHt=zeros(nHt,ntoi,nvox);
allHo=zeros(nHo,ntoi,nvox);

disp('Loading hetero brain data...')
for i=1:nHt
    disp([num2str(i) ' out of ' num2str(nHt)])
    nii=load_nii(cfg.infilesHt{i});
    temp=permute(nii.img,[4 1 2 3]);
    allHt(i,:,:)=zscore(temp(cfg.toi,inmask));
end

disp('Loading homo brain data...')
for i=1:nHo
    disp([num2str(i) ' out of ' num2str(nHo)])
    nii=load_nii(cfg.infilesHo{i});
    temp=permute(nii.img,[4 1 2 3]);
    allHo(i,:,:)=zscore(temp(cfg.toi,inmask));
end

fprintf('\n')


%% Average and save data

% average over voxels
meanHt=mean(allHt,3);
meanHo=mean(allHo,3);

% take 1st principal component over voxels
for i=1:nHt
    temp=squeeze(allHt(i,:,:));
    [~,score]=pca(temp);
    pcaHt(i,:)=score(:,1);
end
for i=1:nHo
    temp=squeeze(allHo(i,:,:));
    [~,score]=pca(temp);
    pcaHo(i,:)=score(:,1);
end

disp('Saving volumes...')

% Check if output directory exists; if not, create it
if ~exist(cfg.outdir,'dir') 
    system(['mkdir -p ' cfg.outdir]);
end

c=regexp(cfg.mask,'[0-9]','match');

filename=[cfg.outdir 'Ht_' cfg.nii(1:end-4) '_brainnetomeROI' c '_meanTC.mat'];
filename=char(strrep(join(filename), ' ', ''));
save(filename,'meanHt');

filename=[cfg.outdir 'Ho_' cfg.nii(1:end-4) '_brainnetomeROI' c '_meanTC.mat'];
filename=char(strrep(join(filename), ' ', ''));
save(filename,'meanHo');

filename=[cfg.outdir 'Ht_' cfg.nii(1:end-4) '_brainnetomeROI' c '_pcaTC.mat'];
filename=char(strrep(join(filename), ' ', ''));
save(filename,'pcaHt');

filename=[cfg.outdir 'Ho_' cfg.nii(1:end-4) '_brainnetomeROI' c '_pcaTC.mat'];
filename=char(strrep(join(filename), ' ', ''));
save(filename,'pcaHo');

disp('Done!')
fprintf('\n')