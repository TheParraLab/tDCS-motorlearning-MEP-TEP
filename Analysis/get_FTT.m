function [ts_all,a_all,ncs_all,r,micro_on,micro_off,noresponse,longpause,nkp] = get_FTT(sid,folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculates FTT performance for a given subject
%   
%   Inputs:
%   sid -       subject ID
%   folder -    directory where behavioral data files are stored
%
%   Outputs:
%   ts_all -    tapping speed for all trials across all task iterations
%   a_all -     accuracy for all trials
%   ncs_all -   number of correct sequences for all trials
%   r -         sensation ratings
%   micro_on -  micro-online performance
%   micro_off - micro-offline performance
%   noresponse -a record of whether there were trials with no response
%   longpause - a record of whether there were trials with long pauses
%   nkp -       total number of keypresses throughout tasks
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(folder) %Folder where raw behavioral data files are stored

sid = num2str(sid); %Subject ID
ntrials = 36;   %Number of trials per section
nsections = 4;  %Number of sections/tasks

ts_all = zeros(ntrials,nsections);      %Tapping speed for all tasks
a_all = zeros(ntrials,nsections);       %Accuracy for all tasks
ncs_all = zeros(ntrials,nsections);     %Number of correct sequences for all tasks
r = zeros(3,1);                         %VAS Sensation ratings
noresponse = zeros(nsections,1);        %Number of trials with no response
longpause = zeros(nsections,1);         %Number of trials with pauses longer than pauselimit
pauselimit = 3;                         %Pause limit in seconds
ises = 1;                               %Session counter: two sessions total: 1 and 2
iseq = 1;                               %Sequence counter: only one for session 1, three different sequences for session 2
nkp = zeros(nsections,1);               %Total number of key presses throughout task

micro_on = zeros(ntrials,nsections); %Micro-online learning for tapping speed
micro_off = zeros(ntrials,nsections); %Micro-offline learning for tapping speed
%to calculate micro-online: subtract speed/interval of first correct sequence from that of the last correct sequence
%to calculate micro-offline: subtract speed/interval of last correct sequence of the previous trial from that of the first correct sequence in the following trial

for m = 1:nsections
    clear log
    session = num2str(ises);
    sn = num2str(iseq);
    filename = strcat('Subject',sid,'_',session,'_S',sn,'.mat');
    load(strcat(folder,filename),'log','sequence') 
    %log is structure array containing a record of keys pressed and the timing of each keystroke, per trial
    %sequence is the 5-element sequence from the task

    if ~exist('log','var')
        ts_all(:,m) = nan;
        a_all(:,m) = nan;
        ncs_all(:,m) = nan;
        disp(strcat(filename," not found"))
    else

        if strcmp(session,'1') && ~isempty(who('-file',filename,'rating'))  %record sensation ratings
            load(strcat(folder,filename),'rating')
            r = rating;
        end
    
        %Keypad keys are mapped to letters on a keyboard: 1=v,2=b,3=n,4=m. Use strcmp to match the keys
        if sequence == ["4";"1";"3";"2";"4"]
            seq = 'mvnbm';
        elseif sequence == ["2";"3";"1";"4";"2"]
            seq = 'bnvmb';
        elseif sequence == ["3";"4";"2";"1";"3"]
            seq = 'nmbvn';
        end
    
        nseq = zeros(ntrials,1);    %number of completed sequences within each trial
        correct = zeros(ntrials,1); %number of correct sequences within each trial
        tt = zeros(ntrials,1);      %array of mean tapping time intervals for all trials
        a = zeros(ntrials,1);       %accuracy
        ttt_last = NaN;             %mean interval of last complete sequence in previous trial
        for i = 1:ntrials
            c = 0; %number of correct presses in one trial
            e = 0; %number of erroneous presses in one trial
    
            logi = log.(strcat('trial',num2str(i))); %load data for the current trial
            if ~isempty(logi)
                key = logi(:,1); %key presses
                time = logi(:,2); %time stamps
                nkp(m) = nkp(m)+length(key); %add to total number of key presses
            else
                key = 0;
                time = 0;
                if i>1  %due to a bug in the task script (now fixed), trial 1 is missing in some samples
                    noresponse(m) = noresponse(m)+1; %count number of trials with no response per section/task
                end
            end
    
            nseq(i) = floor(length(log.(strcat('trial',num2str(i))))/5); %split key presses into sequences of 5 presses
            ttt = 0; %array of tapping time intervals within one trial
            if nseq(i) > 0
                for j = 0:nseq(i)-1
                    cs = 0; %number of correct presses within one sequence
                    for k = 1:5
                        if strcmp(key(j*5+k),seq(k)) %check for correctness at each position within the sequence
                            c = c+1;
                            cs = cs+1;
                        else
                            e = e+1;
                        end
                    end
                    if cs == 5 %check for complete match within sequence
                        t = mean(diff(str2double(time(j*5+1:j*5+5)))); %mean time interval within sequence
                        ttt = [ttt;t]; %speed is counted only for correct sequences
                        correct(i) = correct(i) + 1; %count perfectly correct sequences
                    end
                end
            end
            ic = mod(size(log.(strcat('trial',num2str(i))),1),5);
            if ic > 1 %check for incomplete correct sequence at end of trial
                cs = 0;
                for k = 1:ic
                    if strcmp(key(nseq(i)*5+k),seq(k)) %check for correctness
                        cs = cs+1;
                    end
                end
                if cs == ic %check for complete match within sequence
                    t = mean(diff(str2double(time(nseq(i)*5+1:nseq(i)*5+ic))));
                    ttt = [ttt;t];
                end
            end
            if length(ttt) > 1
                ttt(1) = [];    %remove leading 0
                
                micro_on(i,m) = 1/ttt(end)-1/ttt(1);
                if i>1 && ~isnan(ttt_last)
                    micro_off(i-1,m) = 1/ttt(1)-1/ttt_last;
                elseif i>1 && isnan(ttt_last)
                    micro_off(i-1,m) = 1/ttt(1);
                end
                ttt_last = ttt(end); %interval of last complete sequence in previous trial
            else
                ttt_last = NaN;
            end
            tt(i) = mean(ttt);
            if c > 0
                a(i) = 1 - e/(e+c);
            else
                a(i) = 0;
            end
    
            kpinterval = diff(str2double(time)); %check for long intervals between key presses
            if (sum(kpinterval>pauselimit)||isempty(logi))
                longpause(m) = longpause(m)+1;
            end
        end
        ts = 1./tt;     %tapping speed in presses/s
        ts(~tt) = 0;     %find empty time intervals
    
        ts_all(:,m) = ts;
        a_all(:,m) = a;
        ncs_all(:,m) = correct;
    end
    if ises ~= 1
        iseq = iseq + 1;
    end
    if ises == 1
        ises = ises + 1;
    end
end