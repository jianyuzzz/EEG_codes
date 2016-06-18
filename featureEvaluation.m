%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Data analysis of dataset SPUELER2015
% Feature evaluation and visualization
% 
% Input:    connectivity information (e.g. correlation matrices)
% Output:   connected topoplots

% Author: Jianyu Zhao
% Last revised: 18.06.2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;

inpath = strcat(pwd,'/data_feat/');
filename = 'S01';
load(strcat(inpath,filename));
N = 10;

cov1 = FEAT.conn(:,:,FEAT.labels==1);
cov0 = FEAT.conn(:,:,FEAT.labels==0);

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
[sortedDf,idx] = sort(df(:),'descend');
maxDfs = sortedDf(1:N);
maxIdx = idx(1:N);
[I,J] = ind2sub(size(df),maxIdx);

%plot(sortedDf);

%gscatter(squeeze(FEAT.conn(I(1),J(1),:)),squeeze(FEAT.conn(I(2),J(2),:)),FEAT.labels)


%% Draw the topoplots
% for execution errors
displayStr.chanPairs = [I J];
mcov1 = mean(cov1,3);
mcov0 = mean(cov0,3);
strth = [];
for i=1:size(I)
    for j=1:size(J)
        strth = [strth;mcov1(I(i),J(j))-mcov0(I(i),J(j))];
    end
end
displayStr.connectStrength = strth;
topoplot_connect(displayStr, FEAT.chanlocs);
clear strth;
%{
xp = min(strth);
yp = max(strth);
label = xp:(yp-xp)/10:yp;
label = round((label*100))/100; % round to 2 decimals
st = cell(1,11);
for i=1:11
    st{i} = num2str(label(i));
end
%}

% for no-error situation
%{
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
%}