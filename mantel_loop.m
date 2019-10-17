
% Mantel test between similarity matrices
% out: Nii_corr.nii   : correlation values
%      Nii_pvalue.nii : p values for correlation values

clear
addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI/');
addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/')
addpath('/m/nbe/scratch/braindata/kauppim1/scripts/clusterstats/')


%% Loading data
datapath='/m/nbe/scratch/braindata/jaalho/narrative/';
% braindata; should contain memMaps structure
% load([datapath 'dynISC/assN16_15TRwin_1TRstep_pausesRegressed/memMaps.mat']);

% time windows corresponding to the nine paragraphs
% tws=[
% 21
% 50
% 91
% 109
% 138
% 180
% 217
% 240
% 268
% ];

% mask
% iscmask=load_nii([datapath 'ass_mask.nii']);
iscmask=load_nii([datapath 'mainResult_mask.nii']);
iscmask=iscmask.img;

% behavioral data
behav_data = load([datapath 'mean_simmats_paras_english_first3Words_LSAwordnetMeanBoost_europarlCorpus_n16.mat']);

resultID='english_first3Words_LSAwordnetMeanBoost_meanSimmat_europarlCorpus_pausesRegressed_n16_mainResultMask';

%% Begin loop over paragraphs
Npara=size(behav_data.mean_simmats_paras,3);
for para=1:Npara
    
    disp(['processing para' num2str(para)]);
    % braindata; should contain memMaps structure
    load([datapath 'ISC_para' num2str(para) '_N16/memMaps.mat']);
    % pick the time window and the corresponding behavioral model
%     tw=['timeInt' num2str(tws(para))];
    model=behav_data.mean_simmats_paras(:,:,para);
    
    %Finding all voxels which contain data
    [x,y,z] = ind2sub(size(iscmask),find(iscmask ==1));

    % Calculating all correlation coefficients within the mask
    outvec = []; %Saving correlation values

    %Saving x,y,z coordinates of correlation values
    X = []; 
    Y = [];
    Z = [];

    %Calculate the side length of brain data square matrix
    syms a;
%     temp=double(squeeze(memMaps.cormatMap.win.band0.Session1.cor.(tw).Data.xyzc(x(1),y(1),z(1),:)));
    temp=double(squeeze(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc(x(1),y(1),z(1),:)));
    n = solve(a^2-a-2*length(temp));
    n = double((n(n>0)));

    ids=find(triu(ones(n),1));
    ids_model=find(triu(ones(double(size(model))),1));

    %% loop
    parfor index=1:length(x)
%     for index=1:length(x)

        %Taking those voxels which contains data
%         temp=double(squeeze(memMaps.cormatMap.win.band0.Session1.cor.(tw).Data.xyzc(x(index),y(index),z(index),:)));
        temp=double(squeeze(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc(x(index),y(index),z(index),:)));

        %Modifying temp
        isc_mat=zeros(n);%

        isc_mat(ids)=temp;
        isc_mat=isc_mat+isc_mat'+eye(n);
        final_isc_mat = isc_mat;

        %Correlation with our model
        out=corr(final_isc_mat(ids_model), model(ids_model),'type', 'Spearman');  

        %Saving coefficients
%         outvol_corr(x(index),y(index),z(index)) = out;
        outvec = [outvec ; out];
        X = [X ; x(index)];
        Y = [Y ; y(index)];
        Z = [Z ; z(index)];

    end
    disp('first parfor loop done!');

    %% Mapping all correlation values into right voxels

    outvol_corr = zeros(91,109,91);
    for indx = 1:length(X)
        c_val = outvec(indx); %Correlation value we are calculating
        x_indx = X(indx); %Coordinates of correlation value
        y_indx = Y(indx);
        z_indx = Z(indx);

        outvol_corr(x_indx, y_indx, z_indx) = c_val;
    end

    % Taking 101 correlation values which cover all correlation values

    %Sorting all correlation values
    outvec_S = sort(outvec); 

    % Then we take 100 values that cover the hole range of the
    % distribution. 
%     corrs = prctile(outvec_S,[0:100]); %the 0 percentile is the minimum of all values and 100 percentile is the maximum
    corrs = prctile(outvec_S,[0:100]);

    % Mapping all 101 corrs to right voxels
    outvol_corr_101 = zeros(91,109,91);

    %x,y,z coordinates of 101 correlation values
    X_101 = []; 
    Y_101 = []; 
    Z_101 = []; 
    corrs_101 = [];

    for i = 1:length(corrs)

        [value place] = min(abs(outvec-corrs(i))); % < 0.0001)

        outvol_corr_101(place) = corrs(i);

        X_101 = [X_101 ; X(place)];
        Y_101 = [Y_101 ; Y(place)];
        Z_101 = [Z_101 ; Z(place)];
        corrs_101 = [corrs_101;corrs(i)];
    end


    % Calculating the pvalues for those 101 correlation values

    pvalue_vol_101 = zeros(91,109,91); %Saving pvals as volume
    pvals = []; %Saving pvals as vector

    [x,y,z] = ind2sub(size(iscmask),find(iscmask ==1));
    voxels = find(iscmask); %voxels

    parfor id = 1:length(corrs) %101 

        cor = corrs_101(id);
        x_est = X_101(id); % x coordinate of correlation
        y_est = Y_101(id); % y coordinate of correlation
        z_est = Z_101(id); % z coordinate of correlation

        %Taking brains
        temp_est=double(squeeze(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc(x_est,y_est,z_est,:))); %ISC of chosen correlation coefficient 
        isc_mat=zeros(n); 
        isc_mat(ids)=temp_est;
        isc_mat=isc_mat+isc_mat'+eye(n);
        final_isc_mat = isc_mat;

        iter = 100000; 
        surro = zeros(iter,1); %Estimated correlation values

        for i = 1:iter
            pe = randperm(size(final_isc_mat,1));
            temp_model = final_isc_mat(pe,pe);
            surro(i) = corr(temp_model(ids_model),model(ids_model), 'type', 'spearman');
        end

        [fi xi] = ksdensity(surro, 'function', 'cdf', 'npoints', 200);

        pvalue_left = interp1([-1 xi 1],[0 fi 1], cor); %Calculating pval for chosen correlation value
        pvalue_right = 1 - pvalue_left;
        pvals = [pvals, pvalue_right]; %p value vector

    end
    disp('second parfor loop done!');

    %Mapping those 101 pvals into right place in volume
    for i=1:length(X_101)
        p = pvals(i); %pval we are calculating
        x=X_101(i); %Coordinates of pval
        y=Y_101(i);
        z=Z_101(i);
        pvalue_vol_101(x,y,z) = p;

    end

    % Calculating all pvals based on the estimated CDF
    pvals1 = sort(pvals',1,'descend');
    % For each voxel, mapping its correlation value to the p value using CDF
    outvol_pval= zeros(91,109,91);
    for vox = 1:length(voxels)
        %disp(vox)
        voxel=voxels(vox);
        pv_out = 1-interp1(corrs,pvals1,outvol_corr(voxel), 'pchip');
        outvol_pval(voxel) = pv_out;
    end

    %% Saving unthresholded results
    %Fixing originator:
    disp('saving results');
    filename=[datapath 'Nii_corr_' resultID '_para' num2str(para) '.nii'];
    save_nii(make_nii(outvol_corr), filename);
    nii = bramila_fixOriginator(filename);
    save_nii(nii,filename);
    filename=[datapath 'Nii_pvalue_' resultID '_para' num2str(para) '.nii'];
    save_nii(make_nii(outvol_pval), filename);
    nii = bramila_fixOriginator(filename);
    save_nii(nii,filename);

    %% FDR
    num_pl = numel(find(outvol_pval(voxels) < 0.05));
    q=mafdr(outvol_pval(voxels), 'BHFDR', 'True');
    voxels_fdr = find(q<0.05);
    numel_pvals = numel(voxels_fdr);
    inmask = find(iscmask);
    significant_voxels = inmask(voxels_fdr);

    new_brain = zeros(91,109,91);
    new_brain(significant_voxels)=outvol_corr(significant_voxels);
    filename=[datapath 'Nii_corr_' resultID '_FDR05_para' num2str(para) '.nii'];
    save_nii(make_nii(new_brain),filename);
    nii=bramila_fixOriginator(filename);
    save_nii(nii,filename);

    %% Cluster
    th_l=-1;
    thresholds = [];
    temp_corr=outvol_corr;
    temp_pval=outvol_pval;

    %Setting a threshold of correlation value for every "brain"
    temp_pval = 1-temp_pval;
    th_r = min(temp_corr(find(temp_pval<0.05 & temp_pval>0)));

    thresholds = [thresholds th_r];

    temp_corr(intersect(find(temp_corr<=th_r),find(temp_corr>=0)))=0;
    temp_corr(intersect(find(temp_corr>=th_l),find(temp_corr<=0)))=0;
    mask=clusterit(abs(sign(temp_corr)),1,5,18);
    resultTH=temp_corr.*sign(iscmask);
    outp=getallstats(resultTH,1);
    outn=getallstats(resultTH,-1);
    filename=[datapath 'Nii_corr_' resultID '_clusterCorr_para' num2str(para)];
    csvwrite([filename '_p.csv'],outp);
    csvwrite([filename '_n.csv'],outn);
    filename=[filename '.nii'];
    save_nii(make_nii(resultTH),filename);
    nii=fixOriginator(filename);
    save_nii(nii,filename);
end

