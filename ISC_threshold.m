%%
clear

addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))  
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/'))

paths={
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_allSubs_1stPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_allSubs_2ndPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_allSubs_fullMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_heteros_1stPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_heteros_2ndPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_heteros_fullMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_homos_1stPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_homos_2ndPartMovie'
'/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_homos_fullMovie'
};

%%
for i=1:length(paths)
    % load data
    load([paths{i} '/memMaps.mat'])
    load([paths{i} '/stats/Thband0Session1win0.mat'])
    % fix the name
    memMaps.cormatMap.whole.band0.Session1.cor.Filename=[paths{i} '/stats/corMat_cor_band0_Session1_DWT.bin'];
    % fisher z-transform and the inverse of that (back to the correlation
    % scale)
    cMat=memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc; 
    ISC=tanh(mean(atanh(cMat(:,:,:,:)),4));
    % Th consists 12 threshold values corresponding to
    % {'none_0.05'}    {'pID_0.05'}    {'pN_0.05'}    {'bonf_0.05'}    
    % {'none_0.01'}    {'pID_0.01'}    {'pN_0.01'}    {'bonf_0.01'}    
    % {'none_0.001'}   {'pID_0.001'}   {'pN_0.001'}   {'bonf_0.001'}
    % thresholding (choose the one you like); e.g. 2 -> pID_0.05
    thval=10;
%     Th_info(thval)
    ISCth=zeros(size(ISC,1),size(ISC,2),size(ISC,3));
    ISCth(find(ISC>Th(thval)))=ISC(find(ISC>Th(thval)));

    % save unthresholded
    % filename=[path subjs{i} '_ISC_unthresholded.nii'];
    % save_nii(make_nii(ISC),filename);
    % nii=bramila_fixOriginator(filename);
    % save_nii(nii,filename);

    % save thresholded
    %filename=[path subjs{i} '/ISC_thresholded.nii'];
    filename=[paths{i} '/ISC_thFDR001.nii'];
    save_nii(make_nii(ISCth),filename);
    nii=bramila_fixOriginator(filename);
    save_nii(nii,filename);
end
disp('done!')