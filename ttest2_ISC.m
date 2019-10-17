clear

datapath='/m/nbe/scratch/braindata/jaalho/gaypriest/';
resultID='ttest2_twISC_heteros_vs_homos';
load([datapath 'BL_ref15cons_36pats_30cons_val28pats_10TRwin1TRstep_4mm/memMaps.mat'])
cfg=[];
cfg.infile=memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc;
cfg.mask='/m/nbe/scratch/psykoosi/masks/MNI152_T1_4mm_brain_mask.nii';
cfg.group_id=[zeros(1,15) ones(1,36) 2*ones(1,30) zeros(1,28)];
cfg.p_val_threshold=0.01;
%cfg.doFisherTransform=0;

%%
results=bramila_ttest2_ISC(cfg);

%% Save 3D niftis of raw t-values
filename=[datapath resultID '_tstat.nii'];
save_nii(make_nii(results.raw_tval_map),filename);
nii=fixOriginator(filename,cfg.mask);
save_nii(nii,filename);

%% Save 3D niftis of FDR-thresholded t-values
filename=[datapath resultID '_tstat_FDR' num2str(cfg.p_val_threshold) '.nii'];
newbrain=zeros(size(results.raw_tval_map,1),size(results.raw_tval_map,2),size(results.raw_tval_map,3));
newbrain(find(results.stats.raw_pval_corrected<0.01))=results.raw_tval_map(find(results.stats.raw_pval_corrected<cfg.p_val_threshold));
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

%% Save 3D niftis of TFCE-thresholded t-values
filename=[datapath resultID '_tstat_TFCE' num2str(cfg.p_val_threshold) '.nii'];
newbrain=zeros(size(results.raw_tval_map,1),size(results.raw_tval_map,2),size(results.raw_tval_map,3));
newbrain(find(results.stats.tfce_pval_corrected<0.01))=results.raw_tval_map(find(results.stats.tfce_pval_corrected<cfg.p_val_threshold));
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);

%% Save 3D niftis of cluster-extend-thresholded t-values
filename=[datapath resultID '_tstat_clusterExtend' num2str(cfg.p_val_threshold) '.nii'];
newbrain=zeros(size(results.raw_tval_map,1),size(results.raw_tval_map,2),size(results.raw_tval_map,3));
newbrain(find(results.stats.cluster_pval_corrected<0.01))=results.raw_tval_map(find(results.stats.cluster_pval_corrected<cfg.p_val_threshold));
save_nii(make_nii(newbrain),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);
