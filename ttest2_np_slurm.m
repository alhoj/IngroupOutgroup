clear
delete logs/*
delete jobs/*

cfg=[];

% slurm parameters
cfg.partition='batch';
cfg.mem='64000';
cfg.time='48:00:00';

%% Make jobs -> one job per time window

% time windows of interest
twi=[
3
21
53
77
95
109
125
139
180
201
239
297
308
333
355
384
395
424
433
450
475
488
510
587
597
623
643
678
695];

counter=0;
%
for i=1:length(twi)
    counter=counter+1;
    cfg.toi=twi(i);
    cfg.outdir='/m/nbe/scratch/braindata/jaalho/gaypriest/ttest2/'; % label of the output folder    
    
    cfg.func='ttest2_np';
    cfg.ind=counter; % this is just the job (and log) index
    function_make_scripts_slurm(cfg)

end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');
