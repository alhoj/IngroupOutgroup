% clear
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/');
addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/');
addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/clusterstats/')

datapath='/m/nbe/scratch/braindata/jaalho/gaypriest/';
mask=load_nii([datapath 'group_mask.nii']);
inmask=find(mask.img);

p_val_threshold=0.05;

%%

resultID=['ttest2_np_heteros_vs_homos_twISC_win' num2str(cfg.toi)];

heteros=load([datapath 'twISC/heteros_fullMovie_10TRwin1TRstep/tw' num2str(cfg.toi) '/cormat.mat']);
homos=load([datapath 'twISC/homos_fullMovie_10TRwin1TRstep/tw' num2str(cfg.toi) '/cormat.mat']);

data_heteros=zeros(nchoosek(size(heteros.cormat,2),2),size(mask.img,1),size(mask.img,2),size(mask.img,3));
data_homos=zeros(nchoosek(size(homos.cormat,2),2),size(mask.img,1),size(mask.img,2),size(mask.img,3));

inds_heteros=find(triu(ones(size(heteros.cormat,2)),1)); % indices of the top triangle of the correlation matrix
inds_homos=find(triu(ones(size(homos.cormat,2)),1));

data_heteros(:,inmask)=heteros.cormat(:,inds_heteros)';
data_homos(:,inmask)=homos.cormat(:,inds_homos)';
data_heteros=permute(data_heteros,[2 3 4 1]);
data_homos=permute(data_homos,[2 3 4 1]);

%%
outvol_t=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
outvol_p=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));

for x = 1:size(mask.img,1)
    disp(x)
    temp1=double(squeeze(data_heteros(x,:,:,:)));
    temp2=double(squeeze(data_homos(x,:,:,:)));
    temp1=reshape(temp1,size(data_heteros,2)*size(data_heteros,3),[]);
    temp2=reshape(temp2,size(data_heteros,2)*size(data_heteros,3),[]);
    goodids=find(var(temp1,[],2)>0);
    temp1(find(temp1==-1))=-1+eps;temp1(find(temp1==1))=1-eps;temp1=atanh(temp1);
    temp2(find(temp2==-1))=-1+eps;temp2(find(temp2==1))=1-eps;temp2=atanh(temp2);
    if(length(goodids)>0)
        stats=bramila_ttest2_np([temp1(goodids,:) temp2(goodids,:)],[ones(1,size(temp1,2)) 2*ones(1,size(temp2,2))],5000);
        out_t=zeros(size(temp1,1),1);
        out_t(goodids)=stats.tvals;
        outvol_t(x,:,:)=reshape(out_t,size(data_heteros,2),size(data_heteros,3));
        out_p=ones(size(temp1,1),1);
        out_p(goodids)=stats.pvals(:,1);
        outvol_p(x,:,:)=reshape(out_p,size(data_heteros,2),size(data_heteros,3));
    end
end


% Save uncorrected 3D nifti volume of t-stats and p-vals
%Save t-stat volumes
%     mask=load_nii('/m/nbe/scratch/psykoosi/masks/MNI152_T1_4mm_brain_mask.nii');
filename=[datapath resultID '_tstat.nii'];
save_nii(make_nii(outvol_t),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

%Save p-val volumes
filename=[datapath resultID '_pval.nii'];
save_nii(make_nii(outvol_p),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

%% FDR correction
% load group ISC volume and use it as a mask
% mask=abs(sign(memMaps.resultMap.whole.band0.Session1.cor.Data.xyz));
% mask=load_nii('/m/nbe/scratch/braindata/jaalho/psykoosi/ensipsykoosi/uusi_data/baseline/mask_N136.nii');
inmask = find(mask.img);
% FDR
pvals = 2*min([outvol_p(inmask) 1-outvol_p(inmask)],[],2);
q=mafdr(pvals, 'BHFDR', 'True');
voxels_fdr = find(q<p_val_threshold);

significant_voxels = inmask(voxels_fdr);

newbrain = zeros(size(data_heteros,1),size(data_heteros,2),size(data_heteros,3));
newbrain(significant_voxels)=outvol_t(significant_voxels);

filename=[datapath resultID '_tstat_FDR' num2str(p_val_threshold) '.nii'];
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

%% Cluster correction
msk=clusterit(abs(sign(newbrain)),1,5,18);
resultTH=newbrain.*sign(msk);
outp=getallstats(resultTH,1);
outn=getallstats(resultTH,-1);
% csvwrite([datapath resultID '_clustpos.csv'],outp);
% csvwrite([datapath resultID '_clustneg.csv'],outn);
filename=[datapath resultID '_FDR' num2str(p_val_threshold) '_clusterCorrected.nii'];
save_nii(make_nii(resultTH),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);
