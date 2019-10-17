clear
% delete logs/*
% delete jobs/*

cfg=[];
rois=215:218; % roi numbers from brainnetome atlas
cfg.nii='epi_movie2.nii'; % define nifti file
cfg.toi=[]; % define time points of interest as vector; if empty, take all time points
cfg.outdir='/m/nbe/scratch/braindata/jaalho/gaypriest/BOLDgroupDifference/';
cfg.func='calculate_BOLDgroupDifference(cfg)';

% slurm parameters
cfg.partition='batch';
cfg.mem='64000';
cfg.time='00:29:59';

%% Make jobs

counter=0;

for i=1:length(rois)
    counter=counter+1;
    cfg.mask=['/m/nbe/scratch/braindata/jaalho/gaypriest/brainnetome_rois_separately/' num2str(rois(i)) '.nii']; % define mask
    cfg.ind=counter; % this is just the job (and log) index
    function_make_scripts_slurm(cfg)
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');
