%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Data analysis of dataset SPUELER2015
% Feature extraction
% 
% Input:    struct file containing 3 types of epoched data, no error, 
%           execution error, outcome error (outcome error is the human
%           error, while execution error is the machine error?)
% Output:   struct file containing different feature types and labels
%
% Author: Stefan Ehrlich
% Modificated by Jianyu Zhao
% Last revised: 23.02.2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all

inpath = strcat(pwd,'/data_epo/');% horizontal string concentenation
outpath = strcat(pwd,'/data_feat/');% owd shows the current working directory

% List of filenames of all subjects
% subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};
subjects = {'S01'};
%%
for s=1:length(subjects)
    %%
    filename = subjects{s};% e.g. 'S01'
    %filename = strcat(subjects{s},'.mat');
    load(strcat(inpath,filename))
    
    FEAT.srate = EPO.srate;
    FEAT.chanlocs = EPO.chanlocs;% sampling rate and channel locations 
                                 % remain the same in feature data files
    
    %% Cosntruct label vector
    nTrials_noError = size(EPO.noError,3);% number of no-error trials
    nTrials_exec = size(EPO.exec,3);
    nTrials_out = size(EPO.out,3);
    % the 3rd dim means the number of corresponding epochs
    
    FEAT.labels = [zeros(nTrials_noError,1);zeros(nTrials_exec,1)+1;zeros(nTrials_out,1)+2];
    % 3 columns with different lengths:[0;1;2]
    
    %% Extraction of temporal features
    FEAT.times = EPO.times;% same selected epochs
    FEAT.temp = cat(3,EPO.noError,EPO.exec,EPO.out);% simply contatenate the epoched data 
    
    %% Extraction of frequency-domain features
    %{
    winSize = 128;
    noverlap = 64;
    nfft = 256;% length of fft 
    
    for ch = 1:size(FEAT.temp,1)% for every channel, 28 
        for trial = 1:size(FEAT.temp,3)% for every trial
            [pxx(:,ch,trial) FEAT.f] = pwelch(FEAT.temp(ch,EPO.times>0.2 & EPO.times<0.9,trial),winSize,noverlap,nfft,EPO.srate);
            % estimate the spectral density of temporary features in this
            % case
        end
        fprintf(strcat(num2str(ch),'\n'))
    end
    FEAT.freq = permute(pxx(1:40,:,:),[2 1 3]);% rearrange the dimensions of A to 2-by-1-by-3
    FEAT.f = FEAT.f(1:40);% only preserve the main 40 frequencies
    clear pxx
    
    figure;
    figure(s) % frequency-domain features of the 8th channel
    plot(FEAT.f,mean(FEAT.freq(8,:,FEAT.labels==0),3))
    hold on
    plot(FEAT.f,mean(FEAT.freq(8,:,FEAT.labels==1),3),'g')
    plot(FEAT.f,mean(FEAT.freq(8,:,FEAT.labels==2),3),'r')
    
    %}
    %% Extraction of connectivity features

    for trial = 1:size(FEAT.temp,3)% correlation
        FEAT.conn(:,:,trial) = corr(FEAT.temp(:,EPO.times>0.2 & EPO.times<0.9,trial)');
    end
   


    for trial = 1:size(FEAT.temp,3) % covariance
        %{
        % normalize
        FEAT.conn2(:,:,trial) = cov(FEAT.temp(:,:,trial)')-(trace(cov(FEAT.temp(:,:,trial)')))/28;
        % standarize
        FEAT.conn2(:,:,trial) = FEAT.conn2(:,:,trial)/std(diag(cov(FEAT.temp(:,:,trial)')));
        %FEAT.conn2(:,:,trial) = cov(FEAT.temp(:,:,trial)');
        %}
        covv = cov(FEAT.temp(:,:,trial)');% for each trial, 28 channels, each with 256 moments
        % 256*28, columnwise covariances
        FEAT.conn2(:,:,trial) = (covv - mean(covv(:)))/std(covv(:));

    end   
    
 %{   
    % plot the differences between 
    figure;
    subplot(2,2,1)
    imagesc(mean(FEAT.conn2(:,:,FEAT.labels==1),3)-mean(FEAT.conn2(:,:,FEAT.labels==0),3))
    title('Execution error minus no error')
    cbar%mean of the 3rd dim, i.e. mean of covariance relationship between channels across different trials
    subplot(2,2,2)
    imagesc(mean(FEAT.conn2(:,:,FEAT.labels==2),3)-mean(FEAT.conn2(:,:,FEAT.labels==0),3))
    title('Output error minus no error')
    cbar
    subplot(2,2,3)
    imagesc(mean(FEAT.conn2(:,:,FEAT.labels==2),3)-mean(FEAT.conn2(:,:,FEAT.labels==1),3))
    title('Output error minus execution error')
    cbar
%}
    %%  coherence calculation
    for trial = 1:size(FEAT.temp,3)% coherence
        FEAT.conn3(:,trial) = mscohere(FEAT.temp(:,4,trial)',FEAT.temp(:,12,trial)');
        % 28*256 -> 256*28, columnwise, still not correct though
        % size: 129*28, 129 is the size of normalized freq
        % I can calculate the 28*28 coherence relation, but which freq can
        % i use to build features?
    end
    
    
    
    save(strcat(outpath,filename),'FEAT')
    %clear FEAT
    
    %% Show the temporary waveforms of Fz, Cz and Pz data after exe-error
%{    
    figure; plot(mean(EPO.out(4,:,:),3));title('Fz channel after outcome error');
    figure; plot(mean(EPO.out(12,:,:),3));title('Cz channel after outcome error');
    figure; plot(mean(EPO.out(20,:,:),3));title('Pz channel after outcome error');
%}    
    figure; plot(squeeze(EPO.out(4,:,70)));title('Fz channel after outcome error');
    figure; plot(squeeze(EPO.out(12,:,70)));title('Cz channel after outcome error');
    figure; plot(squeeze(EPO.out(20,:,70)));title('Pz channel after outcome error');

    
    %% coherence plot
    show =  FEAT.conn3(:,1802);
    plot(show);
end
    
    
    
    