clear 
addpath(genpath('/m/nbe/scratch/psykoosi/scripts/NIFTI'));

subs={


};

% cbr=load_nii('/m/nbe/scratch/braindata/jaalho/psykoosi/ensipsykoosi/uusi_data/Cerebrum_mask.nii');

mask=load_nii(subs{1});
disp(1)
for k=2:length(subs)
    temp=load_nii(subs{k});
    mask.img=mask.img.*temp.img;
%     mask.img=mask.img.*temp.img.*cbr.img;
    disp(k)
end

save_nii(mask,'./groupMask.nii');

