clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

cfg.subs = 'homos';
cfg.refs = 'heteros';

cfg.nii='epi_movie.nii';

cfg.mask='/m/nbe/scratch/braindata/jaalho/gaypriest/group_mask.nii'; 
cfg.outdir='/m/nbe/scratch/braindata/jaalho/gaypriest/ISC_homos_heterosRef'; % label of the output folder

%%
cfg=calculate_indepRef_ISC(cfg);
