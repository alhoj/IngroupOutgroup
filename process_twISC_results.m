clear 
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI');
addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/');

% resultIDs={'1stPartMovie_10TRwin1TRstep','2ndPartMovie_10TRwin1TRstep','fullMovie_10TRwin1TRstep'};
resultIDs={'fullMovie_10TRwin1TRstep'};
datapath='/m/nbe/scratch/braindata/jaalho/gaypriest/twISC/';
mask=load_nii('/m/nbe/scratch/braindata/jaalho/gaypriest/group_mask.nii');
inmask=find(mask.img);

nHt=15; % number of heteros
nHo=14; % number of homos
%%
for i=1:length(resultIDs)
    disp(resultIDs{i})
    outlabel=['twISC_' resultIDs{i}];
    ntw=length(dir([datapath 'heteros_' resultIDs{i} '/tw*'])); % number of time windows
    
    dataHt=zeros(length(inmask),ntw); % empty voxels x time windows matrix for group average twISC
    dataHo=zeros(length(inmask),ntw);
    
    tstats=zeros(length(inmask),ntw); % empty voxels x time windows matrix for t-statistics
    
    indsHt=find(triu(ones(nHt),1)); % indices of the top triangle of the correlation matrix
    indsHo=find(triu(ones(nHo),1));
    
    for t=1:ntw
        disp(t);
        tempHt=load([datapath 'heteros_' resultIDs{i} '/tw' num2str(t) '/cormat.mat']);
        tempHo=load([datapath 'homos_' resultIDs{i} '/tw' num2str(t) '/cormat.mat']);
 
        dataHt(:,t)=tanh(mean(atanh(tempHt.cormat(:,indsHt)),2)); % average over ISC pairs after fisher's z-transform
        dataHo(:,t)=tanh(mean(atanh(tempHo.cormat(:,indsHo)),2));
        
        [~,~,~,stats]=ttest2(tempHt.cormat(:,indsHt)',tempHo.cormat(:,indsHo)');
        tstats(:,t)=stats.tstat;       
    end

    fnHt=['/m/nbe/scratch/braindata/jaalho/gaypriest/heteros_' outlabel];
    fnHo=['/m/nbe/scratch/braindata/jaalho/gaypriest/homos_' outlabel];    
    save(fnHt,'dataHt');
    save(fnHo,'dataHo');
    
    % save hetero 4D nifti
    newbrain4D=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),ntw);
    newbrain4D=permute(newbrain4D,[4 1 2 3]);
    newbrain4D(:,inmask)=dataHt';
    newbrain4D=permute(newbrain4D,[2 3 4 1]);
    fnHt=[fnHt '_4D.nii'];
    save_nii(make_nii(newbrain4D),fnHt);
    nii=bramila_fixOriginator(fnHt);
    save_nii(nii,fnHt);
    
    % save homo 4D nifti
    newbrain4D=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),ntw);
    newbrain4D=permute(newbrain4D,[4 1 2 3]);
    newbrain4D(:,inmask)=dataHo';
    newbrain4D=permute(newbrain4D,[2 3 4 1]);
    fnHo=[fnHo '_4D.nii'];
    save_nii(make_nii(newbrain4D),fnHo);
    nii=bramila_fixOriginator(fnHo);
    save_nii(nii,fnHo);
      
    % save hetero vs homo t-test 4D nifti
    newbrain4D=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3),ntw);
    newbrain4D=permute(newbrain4D,[4 1 2 3]);
    newbrain4D(:,inmask)=tstats';
    newbrain4D=permute(newbrain4D,[2 3 4 1]);       
    fn=['/m/nbe/scratch/braindata/jaalho/gaypriest/heteros_vs_homos_' outlabel '_tstats_4D.nii'];
    save_nii(make_nii(newbrain4D),fn);
    nii=bramila_fixOriginator(fn);
    save_nii(nii,fn);
    
end
disp('done!')