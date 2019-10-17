clear 
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/')) 
%% run the nonparametric GLM
load sentiment_scores
regressors={'max_sentiment_scores','min_sentiment_scores'};

cfg=[];
cfg.cdtP = 0.01;
%   cft.cdtR = r value of cluster defining threshold (overrides cdtP)
cfg.alpha = 0.05; % alpha level of the permutation test
cfg.clusterstat = 'maxsum'; % 'maxsize' or 'maxsum'
cfg.NPERM = 5000;
cfg.seed = 0; % seed for the random
cfg.CR=[0 1]; % range to estimate correlations

for r=1:length(regressors)
    for s=1:length(heteroID)
        cfg.infile=['/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Hetero-subjects-rawdata/subject_' num2str(heteroID(s)) '_HT_Prepro/epi_movie.nii'];
        if r==1
            cfg.regressor = max_scores_hetero_regressor_convHRF(:,s); 
        else
            cfg.regressor = min_scores_hetero_regressor_convHRF(:,s);
        end
        cfg.toi = 1:size(cfg.regressor,1); % time points of interest
        cfg=bramila_glm_np(cfg);

        % store the results into a new nifti file       
%         out=cfg.vol(:,:,:,r).*sign(cfg.cmask(:,:,:,r)); % thresholded volume
        out=cfg.vol(:,:,:,r).*sign(cfg.cmask(:,:,:,r)); % unthresholded correlation volume
        filename=['glm/sub' num2str(heteroID(s)) '_bold_vs_' regressors{r} '_corr.nii'];
        save_nii(make_nii(out),filename);
        nii=bramila_fixOriginator(filename);
        save_nii(nii,filename);
    end
end
