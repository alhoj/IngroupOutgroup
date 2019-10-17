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

min_scores_hetero=zeros(str2double(table.time(end)),length(heteroID));
max_scores_hetero=zeros(str2double(table.time(end)),length(heteroID));
min_scores_homo=zeros(str2double(table.time(end)),length(homoID));
max_scores_homo=zeros(str2double(table.time(end)),length(homoID));

%% read the scores into timepoints x subjects matrix
for t=1:str2double(table.time(end))
    disp(t)
    for h=1:length(heteroID)
        try
            min_scores_hetero(t,h)=min(table.negative(table.subject==string(heteroID(h)) & str2double(table.time)==t));
            max_scores_hetero(t,h)=max(table.positive(table.subject==string(heteroID(h)) & str2double(table.time)==t));
        catch
        end
    end
    for h=1:length(homoID)
        
        try
            min_scores_homo(t,h)=min(table.negative(table.subject==string(homoID(h)) & str2double(table.time)==t));
            max_scores_homo(t,h)=max(table.positive(table.subject==string(homoID(h)) & str2double(table.time)==t));
        catch
        end
    end
end
%% average over group
mean_min_scores_hetero=mean(min_scores_hetero,2);
mean_max_scores_hetero=mean(max_scores_hetero,2);
mean_min_scores_homo=mean(min_scores_homo,2);
mean_max_scores_homo=mean(max_scores_homo,2);

%% two-sample t-tests to obtain t-stat timecourses
[~,~,~,stats_min_scores]=ttest2(min_scores_homo',min_scores_hetero');
[~,~,~,stats_max_scores]=ttest2(max_scores_homo',max_scores_hetero');
tstat_min_scores_homo_vs_hetero=stats_min_scores.tstat;
tstat_max_scores_homo_vs_hetero=stats_max_scores.tstat;
tstat_min_scores_hetero_vs_homo=tstat_min_scores_homo_vs_hetero*-1;
tstat_max_scores_hetero_vs_homo=tstat_max_scores_homo_vs_hetero*-1;

%% match the t-stat timecourses with the fMRI volumes to obtain regressors
load tps 
tr=1.7; % repetition time
nvol=714; % number of volumes
trs=tr:tr:nvol*tr; % times of (the endings of) each volume in seconds

data={
max_scores_hetero
min_scores_hetero
max_scores_homo
min_scores_homo
tstat_min_scores_homo_vs_hetero'
tstat_max_scores_homo_vs_hetero'
tstat_min_scores_hetero_vs_homo'
tstat_max_scores_hetero_vs_homo'
max_scores_movieaction
min_scores_movieaction
max_scores_movieother
min_scores_movieother
max_scores_movieprotagonist
min_scores_movieprotagonist
};

labels={
'max_scores_hetero'
'min_scores_hetero'
'max_scores_homo'
'min_scores_homo'
'tstat_min_scores_homo_vs_hetero'
'tstat_max_scores_homo_vs_hetero'
'tstat_min_scores_hetero_vs_homo'
'tstat_max_scores_hetero_vs_homo'
'max_scores_movieaction'
'min_scores_movieaction'
'max_scores_movieother'
'min_scores_movieother'
'max_scores_movieprotagonist'
'min_scores_movieprotagonist'
};
%%
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
save('sentiment_scores_new','sentiment_scores')