%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Data analysis of dataset SPUELER2015
% Epoch extraction
%
% Input:    EEGLAB like data format 
% Output:   struct file containing 3 types of epoched data, no error, 
%           execution error, outcome error
%
% Author: Stefan Ehrlich
% Modificated by Jianyu Zhao
% Last revised: 03.06.2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

inpath = strcat(pwd,'/data_pp/');
outpath = strcat(pwd,'/data_epo/');

% List of filenames of all subjects
% subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};
subjects = {'S01'};

for s = 1:length(subjects)
    
    filename = strcat(subjects{s},'.set');% e.g. S1.set
    out_filename = strcat(subjects{s})% File name for output: S1
    EEG = pop_loadset(strcat(inpath,filename))
    % Load an EEG dataset, named as simple as EEG

    %% copy necessary information
    EPO.srate = EEG.srate;% sampling rate
    EPO.chanlocs = EEG.chanlocs;% channel location information (10-10,10-20)
    
    %% extract trials
    
    [EEG_epo] = pop_epoch(EEG,{'1'},[0 1])
    % convert a continuous EEG dataset to epoched data,
    % by extracting data epochs time locked to specified event types
    EPO.exec = EEG_epo.data;
    [EEG_epo] = pop_epoch(EEG,{'2'},[0 1])% (EEG, events, timelimits)
    EPO.out = EEG_epo.data;% epoched execution error and outcome error data
    
    % recover latencies 
    % record all the moments when errors happend
    execErr = [];
    outErr = [];
    for e = 1:length(EEG.event)
        if EEG.event(e).type == '1'
            execErr = cat(2,execErr,EEG.event(e).latency);
            % {'1'} means there is execution error
            % add the latencies of execution error to the list
        end
        if EEG.event(e).type == '2'
            outErr = cat(2,outErr,EEG.event(e).latency);
            % {'2'}: outcome error
        end
    end
    
    EPO.times = 0:1/EEG.srate:1-1/EEG.srate;
    % a grid from 0 to 1 - 1/EEG.srate; 1/srate means normalized freq.


    % extract noError trials by sliding window
    windowSize = EEG.srate;
    windowShift = EEG.srate;
    errors = cat(2,execErr,outErr);% concatenate the error moments
    c = 1;
    for trial = 1:size(EEG.data,2)/(windowShift)-windowSize
        %## exclude samples near execution and outcome errors
        errortrial = 0;
        for i = 1:length(errors)
            if ((trial-1)*windowShift+windowSize/2 >= errors(i)-EEG.srate) &&...
                ((trial-1)*windowShift+windowSize/2 <= errors(i)+EEG.srate*3)
                errortrial = 1;
            else
                continue;
            end
        end
    
        if errortrial == 0
            EPO.noError(:,:,c) = EEG.data(:,(trial-1)*windowShift+1:(trial-1)*windowShift+windowSize);
            c = c+1;
        end
    end
%%
    figure;
    figure(s)% 3 lines of means of no error, execution error and outcome error epochs
    %plot(EPO.times,mean(EPO.noError(8,:,:),3))
    %hold on
    %plot(EPO.times,mean(EPO.exec(8,:,:),3),'g')
    %plot(EPO.times,mean(EPO.out(8,:,:),3),'r')
    plot(EPO.times,mean(EPO.out(8,:,:),3),'r',EPO.times,...
        mean(EPO.noError(8,:,:),3),'b',EPO.times,mean(EPO.exec(8,:,:),3),...
        'g','LineWidth',0.75);
    legend('outcome error','no error','execution error');
    
    xlabel('Time(s)');
    ylabel('Amplitude(\muV)');
    title('Epoched Data');

    
    %eeglab2loreta(EPO.chanlocs, mean(EPO.exec(:,:,:), 3), 'exporterp', 'on');
   
    
save(strcat(outpath,out_filename),'EPO')

end
