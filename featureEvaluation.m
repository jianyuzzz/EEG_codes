%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Data analysis of dataset SPUELER2015
% Feature evaluation and visualization
% 
% Input:    connectivity information (e.g. covariance matrices)
% Output:   connected topoplots

% Author: Jianyu Zhao
% Last revised: 16.06.2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;

inpath = strcat(pwd,'/data_feat/');
filename = 'S01';
load(strcat(inpath,filename));
N = 10;

% there's no need to calculate the mean value, 
% we need the whole distribution
%{ 
cov1 = mean(FEAT.conn2(:,:,FEAT.labels==1),3); % mean of cov of exe err
cov2 = mean(FEAT.conn2(:,:,FEAT.labels==2),3); % output error
cov0 = mean(FEAT.conn2(:,:,FEAT.labels==0),3); % no error
%}

cov1 = FEAT.conn2(:,:,FEAT.labels==1);
cov0 = FEAT.conn2(:,:,FEAT.labels==0);

n_ch = size(FEAT.temp,1);
[d dp df dfp] = deal(zeros(n_ch,n_ch));

%% Discriminating ability between exe error and no error
for i=1:n_ch
    for j=i+1:n_ch
        [d(i,j) dp(i,j) df(i,j) dfp(i,j)] = ...
            feature_eval_nParametric(cov1(i,j,:),cov0(i,j,:));
    end
end

%% Choose the best N features with the highest df
[sortedDf,idx] = sort(dfp(:),'descend');
maxDfs = sortedDf(1:N);
maxIdx = idx(1:N);
[I,J] = ind2sub(size(df),maxIdx);

%% Draw the topoplots
% for execution errors
displayStr.chanPairs = [I J];
mcov1 = mean(cov1,3);
strth = [];
for i=1:size(I)
    for j=1:size(J)
        strth = [strth;mcov1(I(i),J(j))];
    end
end
displayStr.connectStrength = strth;
topoplot_connect(displayStr, FEAT.chanlocs);
clear strth;
% for no-error situation
mcov0 = mean(cov0,3);
strth = [];
for i=1:size(I)
    for j=1:size(J)
        strth = [strth;mcov0(I(i),J(j))];
    end
end
displayStr.connectStrength = strth;
figure;
topoplot_connect(displayStr, FEAT.chanlocs);