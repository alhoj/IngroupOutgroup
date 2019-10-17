clear

addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/')

load subIDs

mask=load_nii('/m/nbe/scratch/braindata/jaalho/gaypriest/group_mask.nii');
inmask=find(mask.img);
inpath='/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/';
%%
temp=load_nii([inpath 'Hetero-subjects-rawdata/subject_' num2str(heteroID(1)) '_HT_Prepro/epi_movie.nii']);
ntps=size(temp.img,4);

group_data_heteros=zeros(length(heteroID),ntps,length(inmask));
for i=1:length(heteroID)
    disp(i)
    infile=load_nii([inpath 'Hetero-subjects-rawdata/subject_' num2str(heteroID(i)) '_HT_Prepro/epi_movie.nii']);
    data=permute(infile.img,[4 1 2 3]);
    group_data_heteros(i,:,:)=zscore(data(:,inmask));
end

group_data_homos=zeros(length(homoID),ntps,length(inmask));
for i=1:length(homoID)
    disp(i)
    infile=load_nii([inpath 'Homo-subjects-rawdata/subject_' num2str(homoID(i)) '_HO_Prepro/epi_movie.nii']);
    data=permute(infile.img,[4 1 2 3]);
    group_data_homos(i,:,:)=zscore(data(:,inmask));
end

%%
group_data_heteros_mean=squeeze(mean(group_data_heteros,1));
group_data_homos_mean=squeeze(mean(group_data_homos,1));

% save hetero mean bold
filename='heteros_bold_4D.nii';
newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),ntps);
newbrain=permute(newbrain,[4 1 2 3]);
newbrain(:,inmask)=group_data_heteros_mean;
newbrain=permute(newbrain,[2 3 4 1]);
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

% save homo mean bold
filename='homos_bold_4D.nii';
newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),ntps);
newbrain=permute(newbrain,[4 1 2 3]);
newbrain(:,inmask)=group_data_homos_mean;
newbrain=permute(newbrain,[2 3 4 1]);
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

% two-sample t-test
[~,~,~,stats]=ttest2(group_data_homos,group_data_heteros);
tstat_vol=zeros(ntps,size(mask.img,1),size(mask.img,2),size(mask.img,3));
tstat_vol(:,inmask)=squeeze(stats.tstat);
tstat_vol=permute(tstat_vol,[2 3 4 1]);

% save tstat volume
filename='ttest2_homos_vs_heteros_tstats.nii';
save_nii(make_nii(tstat_vol),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);
disp('done!');