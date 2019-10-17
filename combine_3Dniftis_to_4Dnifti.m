clear

addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/'))  
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/'))

heteroIDs={
'sub2'
'sub3'
'sub6'
'sub7'
'sub8'
'sub9'
'sub10'
'sub11'
'sub12'
'sub13'
'sub14'
'sub15'
'sub16'
'sub18'
'sub19'
};

homoIDs={
'sub4'
'sub5'
'sub17'
'sub20'
'sub21'
'sub22'
'sub23'
'sub24'
'sub25'
'sub26'
'sub27'
'sub28'
'sub29'
'sub30'
};

N_hetero=length(heteroIDs);
N_homo=length(homoIDs);
N_all=N_hetero+N_homo;
path='/m/nbe/scratch/braindata/jaalho/gaypriest/';
groups={'all','heteros','homos'};
scores={'face1stPart'};%,'maxSentimentMovie'};

%%
mask=load_nii([path  'group_mask.nii']);
inmask=find(mask.img);
for g=1:length(groups)
    for s=1:length(scores)
        if isequal(groups{g},'all')
            data4D=zeros(N_all,size(mask.img,1),size(mask.img,2),size(mask.img,3));
            for i=1:N_all
                if i <= N_hetero
                    data=load_nii([path '/glm/' heteroIDs{i} '_bold_vs_' scores{s} '_corr.nii']);
                    data4D(i,inmask)=data.img(inmask);
                else
                    data=load_nii([path 'glm/' homoIDs{i-N_hetero} '_bold_vs_' scores{s} '_corr.nii']);
                    data4D(i,inmask)=data.img(inmask);
                end
            end
            data4D=permute(data4D,[2 3 4 1]);
            filename=[path 'glm/all_bold_vs_' scores{s} '_4D.nii'];
            save_nii(make_nii(data4D),filename);
            nii=bramila_fixOriginator(filename);
            save_nii(nii,filename);
        elseif isequal(groups{g},'heteros')
            data4D=zeros(N_hetero,size(mask.img,1),size(mask.img,2),size(mask.img,3));
            for i=1:N_hetero
                data=load_nii([path '/glm/' heteroIDs{i} '_bold_vs_' scores{s} '_corr.nii']);
                data4D(i,inmask)=data.img(inmask);
            end
            data4D=permute(data4D,[2 3 4 1]);
            filename=[path 'glm/heteros_bold_vs_' scores{s} '_4D.nii'];
            save_nii(make_nii(data4D),filename);
            nii=bramila_fixOriginator(filename);
            save_nii(nii,filename);
        elseif isequal(groups{g},'homos')
            data4D=zeros(N_homo,size(mask.img,1),size(mask.img,2),size(mask.img,3));
            for i=1:N_homo
                data=load_nii([path 'glm/' homoIDs{i} '_bold_vs_' scores{s} '_corr.nii']);
                data4D(i,inmask)=data.img(inmask);
            end
            data4D=permute(data4D,[2 3 4 1]);
            filename=[path '/glm/homos_bold_vs_' scores{s} '_4D.nii'];
            save_nii(make_nii(data4D),filename);
            nii=bramila_fixOriginator(filename);
            save_nii(nii,filename);
        end       
    end
end
disp('done!')