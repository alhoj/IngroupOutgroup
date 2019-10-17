clear
% delete logs/*
% delete jobs/*

cfg=[];
cfg.subs='heteros'; % 'heteros', 'homos', or 'both'
cfg.mask='/m/nbe/scratch/braindata/jaalho/gaypriest/group_mask.nii'; % define mask
cfg.nii='epi_movie2.nii'; % define nifti file

% slurm parameters
cfg.partition='batch';
cfg.mem='64000';
cfg.time='00:09:59';

%% Make jobs -> one job per time window
if isequal(cfg.nii,'epi_movie.nii')
    ntps=714; % time points in the original niftis
elseif isequal(cfg.nii,'epi_movie1.nii') || isequal(cfg.nii,'epi_movie2.nii')
    ntps=357; % time points in the original niftis
end
twl=10; % time window length
step=1; % step between consecutive time windows
counter=0;

for i=1:step:ntps
    counter=counter+1;
    tw=i:(twl+i-1);
    cfg.toi=tw; % time points of interest as vector
    if isequal(cfg.nii,'epi_movie.nii')
        cfg.outdir=['/m/nbe/scratch/braindata/jaalho/gaypriest/twISC/' cfg.subs '_fullMovie_' num2str(twl) 'TRwin' num2str(step) 'TRstep'  '/tw' num2str(counter)]; % label of the output folder    
    elseif isequal(cfg.nii,'epi_movie1.nii')
        cfg.outdir=['/m/nbe/scratch/braindata/jaalho/gaypriest/twISC/' cfg.subs '_1stPartMovie_' num2str(twl) 'TRwin' num2str(step) 'TRstep'  '/tw' num2str(counter)]; % label of the output folder
    elseif isequal(cfg.nii,'epi_movie2.nii')
        cfg.outdir=['/m/nbe/scratch/braindata/jaalho/gaypriest/twISC/' cfg.subs '_2ndPartMovie_' num2str(twl) 'TRwin' num2str(step) 'TRstep'  '/tw' num2str(counter)]; % label of the output folder
    end
    cfg.func='calculate_twISC(cfg)';
    cfg.ind=counter; % this is just the job (and log) index
    function_make_scripts_slurm(cfg)
    if tw(end)>=ntps
        break
    end
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');
