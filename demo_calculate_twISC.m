% clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define subjects of interest here
cfg.subs = {
'EPJO203'
'EPJO204'
'EPJO205'
'EPJO206'
'EPJO207'
'EPKK301'
'EPKK302'
'EPPE103'
'EPPE105'
'EPPE106'
'EPPE108'
'EPPE109'
'EPPE112'
};

% Define reference group here
cfg.refs = {
'EPJO203'
'EPJO204'
'EPJO205'
'EPJO206'
'EPJO207'
'EPKK301'
'EPKK302'
'EPPE103'
'EPPE105'
'EPPE106'
'EPPE108'
'EPPE109'
'EPPE112'
};

cfg.condSubs='BL'; % which data to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % which data to use for reference subjects
cfg.res='8mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=[]; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_4mm_brain_mask.nii"
cfg.useMeanOverRefs=0; % take average over the reference group and use that to calculate ISC; options 1->yes or 0->no 

%%
notps=245; % time points in the original niftis
twl=10; % time window length
step=1; % step between consecutive time windows
for i=1:step:notps
    tw=i:twl+i-1;
    if tw(end)>=notps
        break
    end
    cfg.toi=tw; % time points of interest as vector
    cfg.outdir=['test/tw' num2str(i)]; % label of the output folder
    cfg=calculate_indepRef_twISC(cfg);
end

