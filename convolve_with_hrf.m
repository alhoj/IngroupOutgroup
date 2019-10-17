clear; close all;
addpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/')

load('ratings.mat')

hrf=bramila_hrf(1.8); % TR as input
ratings_convHRF=zeros(size(ratings,1),size(ratings,2));
for r=1:size(ratings,2)
    temp=conv(ratings_zscored(:,r),hrf);
    ratings_convHRF(:,r)=temp(1:size(ratings,1));
end
ratings_convHRF_zscored=zscore(ratings_convHRF);