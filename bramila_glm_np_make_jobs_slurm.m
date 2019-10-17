clear

delete logs/*
delete jobs/*

load subIDs
% regressors={'hetero_vs_homo_maxSentimentScores_twISC','homo_vs_hetero_maxSentimentScores_twISC','hetero_vs_homo_minSentimentScores_twISC','homo_vs_hetero_minSentimentScores_twISC'};
% regressors={'maxSentimentScores','minSentimentScores'};
regressors={'face1stPart'};

cfg=[];
cfg.ind=1;
cfg.mask='/m/nbe/scratch/braindata/jaalho/gaypriest/group_mask.nii';
% cfg.toi=358:714;
cfg.func = 'bramila_glm_np_slurm_script';

% slurm parameters
cfg.partition='batch';
cfg.mem='48000';
cfg.time='00:09:59';

inpath='/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/';
% inpath='/m/nbe/scratch/braindata/jaalho/gaypriest/';

% this is for the time-windowed ISC vs sentiment score group-level glm
% for r=1:length(regressors)
%     cfg.infile=[inpath 'heteros_vs_homos_twISC_fullMovie_10TRwin1TRstep_tstats_4D.nii'];
%     cfg.regressor=regressors{r};
%     function_make_scripts_slurm(cfg)
%     cfg.ind=cfg.ind+1;
% end
% for r=1:length(regressors)
%     cfg.infile=[inpath 'homos_vs_heteros_twISC_fullMovie_10TRwin1TRstep_tstats_4D.nii'];
%     cfg.regressor=regressors{r};
%     function_make_scripts_slurm(cfg)
%     cfg.ind=cfg.ind+1;
% end

% for r=1:length(regressors)
%     cfg.infile=[inpath 'ttest2_bold_homos_vs_heteros_tstats.nii'];
%     cfg.regressor=regressors{r};
%     function_make_scripts_slurm(cfg)
%     cfg.ind=cfg.ind+1;
% end
% for r=1:length(regressors)
%     cfg.infile=[inpath 'ttest2_bold_heteros_vs_homos_tstats.nii'];
%     cfg.regressor=regressors{r};
%     function_make_scripts_slurm(cfg)
%     cfg.ind=cfg.ind+1;
% end

% this is for running individual correlation maps; i.e. first level glm
for r=1:length(regressors)
    for s=1:length(heteroID)
        cfg.subi=s;
        cfg.infile=[inpath 'Hetero-subjects-rawdata/subject_' num2str(heteroID(s)) '_HT_Prepro/epi_movie.nii'];
        cfg.regressor=regressors{r};
        function_make_scripts_slurm(cfg)
        cfg.ind=cfg.ind+1;
    end
end
for r=1:length(regressors)
    for s=1:length(homoID)
        cfg.subi=s;
        cfg.infile=[inpath 'Homo-subjects-rawdata/subject_' num2str(homoID(s)) '_HO_Prepro/epi_movie.nii'];
        cfg.regressor=regressors{r};
        function_make_scripts_slurm(cfg)
        cfg.ind=cfg.ind+1;
    end
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');
