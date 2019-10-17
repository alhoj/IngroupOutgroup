function cfg=calculate_indepRef_ISC(cfg)
% Calculates individual ISC maps with respect to a reference group
% 
% Usage:
%   cfg = calculate_indepRef_ISC(cfg);
%

%%
if ~ismember(cfg.subs,{'heteros','homos','both'})
    error('cfg.subs should be either ''heteros'', ''homos'' or ''both''!')
end
if ~ismember(cfg.nii,{'epi_movie.nii','epi_movie1.nii','epi_movie2.nii'})
    error('cfg.nii should be either ''epi_movie.nii'', ''epi_movie1.nii'' or ''epi_movie2.nii''!')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['Could not find mask: ' cfg.mask])
end

disp(['Calculating ISC for ' cfg.subs ' with respect to ' cfg.refs])

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
        cfg.subs=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii);
        subs=heteros;
    case 'homos'
        cfg.subs=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii);
        subs=homos;
    case 'both'
        cfg.subs=[strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii) strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii)];
        subs=both;
end

switch cfg.refs
    case 'heteros'
        cfg.refs=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii);
    case 'homos'
        cfg.refs=strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii);
    case 'both'
        cfg.refs=[strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/', heteros, '_HT_Prepro/', cfg.nii) strcat('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Homo-subjects-rawdata/', homos, '_HO_Prepro/', cfg.nii)];
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
test_nii=load_nii(cfg.subs{1});

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects
ntps=size(temp(:,inmask),1);

nvox=size(temp(:,inmask),2);
nsub=length(cfg.subs);
nref=length(cfg.refs);

% Preallocate data matrices with the dimensions of time points of interest x subjects x voxels 
allsubs=zeros(ntps,nsub,nvox);
allrefs=zeros(ntps,nref,nvox);

disp('Loading brain data of subjects of interest...')
for i=1:nsub
    disp([num2str(i) ' out of ' num2str(nsub)])
    nii=load_nii(cfg.subs{i});
    temp=permute(nii.img,[4 1 2 3]);
    allsubs(:,i,:)=zscore(temp(:,inmask));
%     allsubs(:,i,:)=temp(:,inmask);
end
fprintf('\n')

disp('Loading brain data of reference subjects...')
for i=1:nref
    disp([num2str(i) ' out of ' num2str(nref)])
    nii=load_nii(cfg.refs{i});
    temp=permute(nii.img,[4 1 2 3]);
    allrefs(:,i,:)=zscore(temp(:,inmask));
%     allrefs(:,i,:)=temp(:,inmask);
end

fprintf('\n')


%% Calculate ISC
% Preallocate the correlation matrix
cors=zeros(nvox,nsub,nref);
    
disp('Calculating ISC...')
for voxi=1:nvox
    if mod(voxi,1000)==0 % Show the status every 1000 voxels
        disp([num2str(voxi) '/' num2str(nvox) ' voxels'])
    end
       
    % Calculate the correlation over time of all subjects-of-interest to all reference subjects
    cors(voxi,:,:)=corr(squeeze(allsubs(:,:,voxi)),squeeze(allrefs(:,:,voxi)));
  
end

fprintf('\n')

%% Average to obtain indivual ISC maps

% Take the average across all reference subjects (first Fisher's
% z-tranform, then averaging, then back to correlation scale)
% avg_cors=squeeze(tanh(mean(atanh(cors),3)));
% or without the Fisher's z-transform
%     avg_cors=squeeze(mean(cors,3));

% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

% Check if output directory exists; if not, create it
dirname=cfg.outdir;
if ~exist(dirname,'dir') 
    system(['mkdir -p ' dirname]);
end

for i=1:nsub
    disp(['Subject ' subs{i} ' - ' num2str(i) ' out of ' num2str(nsub)])
    % Take the average across all reference subjects (first Fisher's
    % z-tranform, then averaging, then back to correlation scale)
    avg_cors=squeeze(tanh(mean(atanh(cors(:,i,:)),3)));
    newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
    newbrain(inmask)=avg_cors;
    
    filename=[dirname '/' subs{i} '.nii'];
    save_nii(make_nii(newbrain),filename);
    nii=fixOriginator(filename,mask);
    save_nii(nii,filename);
end

cfg.cors=cors;
disp('Done!')
fprintf('\n')