% this script you don't need to touch unless you want to change some parameters,
% such as the significance thresholds, number of permutations etc.

addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/'))

% run the nonparametric GLM
disp('running GLM');

cfg.cdtP = 0.01;
%   cft.cdtR = r value of cluster defining threshold (overrides cdtP)
cfg.alpha = 0.05; % alpha level of the permutation test
cfg.clusterstat = 'maxsum'; % 'maxsize' or 'maxsum'
cfg.NPERM = 0;
cfg.seed = 0; % seed for the random
cfg.CR=[0 1]; % range to estimate correlations

% load sentiment_scores
load sentiment_scores_new
load subIDs
regressorID=cfg.regressor;
switch regressorID
    case {'face','face1stPart','face2ndPart'}
        load face_regressor
        cfg.regressor = regressor;
    case 'maxSentimentMovie'   
        cfg.regressor = sentiment_scores.max_scores_movie.regressor_convHRF;
    case 'minSentimentMovie'   
        cfg.regressor = sentiment_scores.min_scores_movie.regressor_convHRF;
    case 'meanSentimentMovie'
        cfg.regressor = sentiment_scores.mean_scores_movie.regressor_convHRF;
    case 'maxSentimentScores'
        if contains(cfg.infile,'Hetero')
%             cfg.regressor = max_scores_hetero_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.max_scores_hetero.regressor_convHRF(:,cfg.subi);
        elseif contains(cfg.infile,'Homo')
%             cfg.regressor = max_scores_homo_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.max_scores_homo.regressor_convHRF(:,cfg.subi);
        end
    case 'minSentimentScores'
        if contains(cfg.infile,'Hetero')
%             cfg.regressor = min_scores_hetero_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.min_scores_hetero.regressor_convHRF(:,cfg.subi);
        elseif contains(cfg.infile,'Homo')
%             cfg.regressor = min_scores_homo_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.min_scores_homo.regressor_convHRF(:,cfg.subi);
        end
    case 'meanSentimentScores'
        if contains(cfg.infile,'Hetero')
%             cfg.regressor = min_scores_hetero_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.mean_scores_hetero.regressor_convHRF(:,cfg.subi);
        elseif contains(cfg.infile,'Homo')
%             cfg.regressor = min_scores_homo_regressor_convHRF(:,cfg.subi);
            cfg.regressor = sentiment_scores.mean_scores_homo.regressor_convHRF(:,cfg.subi);
        end
    case 'hetero_vs_homo_maxSentimentScores'
%         cfg.regressor = tstat_max_scores_hetero_vs_homo_regressor_convHRF;
        cfg.regressor = sentiment_scores.tstat_max_scores_hetero_vs_homo.regressor_convHRF;
    case 'hetero_vs_homo_minSentimentScores'
%         cfg.regressor = tstat_min_scores_hetero_vs_homo_regressor_convHRF;
        cfg.regressor = sentiment_scores.tstat_min_scores_hetero_vs_homo.regressor_convHRF;
    case 'homo_vs_hetero_maxSentimentScores'
%         cfg.regressor = tstat_max_scores_homo_vs_hetero_regressor_convHRF;
        cfg.regressor = sentiment_scores.tstat_max_scores_homo_vs_hetero.regressor_convHRF;
    case 'homo_vs_hetero_minSentimentScores'
%         cfg.regressor = tstat_min_scores_homo_vs_hetero_regressor_convHRF;
        cfg.regressor = sentiment_scores.tstat_min_scores_homo_vs_hetero.regressor_convHRF;
    case 'hetero_vs_homo_maxSentimentScores_twISC'
%         cfg.regressor = tstat_max_scores_hetero_vs_homo_regressor_convHRF(5:end-5);
        cfg.regressor = sentiment_scores.tstat_max_scores_hetero_vs_homo.regressor_convHRF(5:end-5);
    case 'hetero_vs_homo_minSentimentScores_twISC'
%         cfg.regressor = tstat_min_scores_hetero_vs_homo_regressor_convHRF(5:end-5);
        cfg.regressor = sentiment_scores.tstat_min_scores_hetero_vs_homo.regressor_convHRF(5:end-5);
    case 'homo_vs_hetero_maxSentimentScores_twISC'
%         cfg.regressor = tstat_max_scores_homo_vs_hetero_regressor_convHRF(5:end-5);
        cfg.regressor = sentiment_scores.tstat_max_scores_homo_vs_hetero.regressor_convHRF(5:end-5);
    case 'homo_vs_hetero_minSentimentScores_twISC'
%         cfg.regressor = tstat_min_scores_homo_vs_hetero_regressor_convHRF(5:end-5);
        cfg.regressor = sentiment_scores.tstat_min_scores_homo_vs_hetero.regressor_convHRF(5:end-5);
end

if isempty(cfg.toi)
    cfg.toi = 1:size(cfg.regressor,1);
end

cfg=bramila_glm_np(cfg);

out=cfg.vol; % unthresholded correlation volume
% out=cfg.vol(:,:,:).*sign(cfg.cmask(:,:,:)); % thresholded volume

if contains(cfg.infile,'Hetero')
    filename=['glm/sub' num2str(heteroID(cfg.subi)) '_bold_vs_' regressorID '_corr.nii']; % unthresholded correlation volume
elseif contains(cfg.infile,'Homo')
    filename=['glm/sub' num2str(homoID(cfg.subi)) '_bold_vs_' regressorID '_corr.nii']; % unthresholded correlation volume
elseif contains(cfg.infile,'bold_homos_vs_heteros')
    filename=['glm/bold_homos_vs_heteros_tstats_VS_' regressorID '_corr_cdt01_perm05.nii'];
elseif contains(cfg.infile,'bold_heteros_vs_homos')
    filename=['glm/bold_heteros_vs_homos_tstats_VS_' regressorID '_corr_cdt01_perm05.nii'];
elseif contains(cfg.infile,'heteros_vs_homos_twISC')
    filename=['glm/twISC_heteros_vs_homos_tstats_VS_' regressorID '_corr_cdt01_perm05.nii'];
elseif contains(cfg.infile,'homos_vs_heteros_twISC')
    filename=['glm/twISC_homos_vs_heteros_tstats_VS_' regressorID '_corr_cdt01_perm05.nii'];
end


save_nii(make_nii(out),filename);
nii=bramila_fixOriginator(filename);
save_nii(nii,filename);
disp('done!');


