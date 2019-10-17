% script for reading and processing data from table table 
% consisting association words and their valence (sentiment analysis score) 
% produced by homo- and heterosexual subjects

clear

addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/')

table=readtable('/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/Sentiment_Analysis/Polarity_table_senti_strength_normalized.csv');
%
heteroID=unique(table.subject(str2double(table.HO)==0)); % codes of heterosexual subejcts
homoID=unique(table.subject(str2double(table.HO)==1)); % codes of homosexual subjects
allID=[heteroID; homoID];

min_scores_all=zeros(str2double(table.time(end)),length(allID));
max_scores_all=zeros(str2double(table.time(end)),length(allID));
%% read the scores into timepoints x subjects matrix
for t=1:str2double(table.time(end))
    disp(t)
    for h=1:length(allID)
        try
            min_scores_all(t,h)=min(table.negative(table.subject==string(allID(h)) & str2double(table.time)==t));
            max_scores_all(t,h)=max(table.negative(table.subject==string(allID(h)) & str2double(table.time)==t));
        catch
        end
    end
end
%%
disp('doing perms')
for iter=1:990
    disp(iter)

    permInds=randperm(length(allID));
    min_scores_group1=min_scores_all(:,permInds(1:length(heteroID)));
    max_scores_group1=max_scores_all(:,permInds(1:length(heteroID)));
    min_scores_group2=min_scores_all(:,permInds(length(heteroID)+1:end));
    max_scores_group2=max_scores_all(:,permInds(length(heteroID)+1:end));
    
    % two-sample t-tests to obtain t-stat timecourses
    
    [~,~,~,stats_min_scores]=ttest2(min_scores_group2',min_scores_group1');
    [~,~,~,stats_max_scores]=ttest2(max_scores_group2',max_scores_group1');
    tstat_min_scores_group2_vs_group1=stats_min_scores.tstat;
    tstat_max_scores_group2_vs_group1=stats_max_scores.tstat;
    tstat_min_scores_group1_vs_group2=tstat_min_scores_group2_vs_group1*-1;
    tstat_max_scores_group1_vs_group2=tstat_max_scores_group2_vs_group1*-1;

    % match the t-stat timecourses with the fMRI volumes to obtain regressors
    load tps 
    tr=1.7; % repetition time
    nvol=714; % number of volumes
    trs=tr:tr:nvol*tr; % times of (the endings of) each volume in seconds

    data={
    max_scores_group1
    min_scores_group1
    max_scores_group2
    min_scores_group2
    tstat_min_scores_group2_vs_group1'
    tstat_max_scores_group2_vs_group1'
    tstat_min_scores_group1_vs_group2'
    tstat_max_scores_group1_vs_group2'
    };

    labels={
    'max_scores_group1'
    'min_scores_group1'
    'max_scores_group2'
    'min_scores_group2'
    'tstat_min_scores_group2_vs_group1'
    'tstat_max_scores_group2_vs_group1'
    'tstat_min_scores_group1_vs_group2'
    'tstat_max_scores_group1_vs_group2'
    };
    %
    for j=1:length(labels)
        % let's flip the negative scores into positive
        if contains(labels{j},'min')
            data{j}=data{j}*-1;
        end

        sentiment_scores.(labels{j}).original=data{j};

        regressor=zeros(nvol,size(data{j},2));

        for i=1:length(tps)
            if i==1
                rep=find(trs<tps(i));            
            else
                rep=find(trs>=tps(i-1) & trs<tps(i));                       
            end
            regressor(rep,:)=repmat(data{j}(i,:),length(rep),1);
        end

        % convolve the regressors with hrf function    
        hrf=bramila_hrf(tr);
        regressor_convHRF=zeros(nvol,size(data{j},2));
        for i=1:size(data{j},2)
            temp=conv(regressor(:,i),hrf);
            regressor_convHRF(:,i)=temp(1:nvol);
        end

        sentiment_scores.(labels{j}).regressor=regressor;
        sentiment_scores.(labels{j}).regressor_convHRF=regressor_convHRF;

    end
    save(['sentiment_scores_perm/sentiment_scores_perm' num2str(iter)],'sentiment_scores')
end