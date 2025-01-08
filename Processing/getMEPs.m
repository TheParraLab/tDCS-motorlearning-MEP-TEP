function [meps,traces,ttrace,fs] = getMEPs(sn,isession,ihand,TMSdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Extracts MEPs for a given subject, session number, and hand
%
%   Inputs:
%   sn -        Subject ID number, e.g. 33
%   isession -  Session index (pre/post = 1/2)
%   ihand -     Hand index (L/R = 1/2)
%   TMSdir -    Directory where raw EMG recordings are stored
%
%   Outputs:
%   meps -      Array of MEP amplitudes
%   traces -    Aligned EMG traces
%   ttrace -    Timestamps corresponding to EMG traces
%   fs -        Sampling rate
%
%   Requires Libeep MATLAB importer (tested with 3.3.177)
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Processing parameters
before = 0.5;       %seconds before trigger
after = 0.5;        %seconds after trigger
fc = 1;             %cutoff frequency for high pass filter
tend = 0.3;         %seconds at the end to sample for noise
f = 58:0.1:62;      %frequency range to test for harmonics
K = 6;              %harmonic analysis order
tmepmin = 0.019;    %MEP start time
tmepmax = 0.05;     %MEP end time

%Load the table containing TMS trial numbers
trialsfile = 'TMStrials.csv';
trials = readtable(strcat(TMSdir,'\',trialsfile));
sid = strcat('A',num2str(sn)); %This was hardcoded for this study. Change accordingly if necessary
itable = find(strcmp(trials.Subject,sid)); %row number matching current subject
%Column labels used in the table
triggerlabels1 = ["LHPreFirst","RHPreFirst";"LHPostFirst","RHPostFirst"];
triggerlabels2 = ["LHPreLast","RHPreLast";"LHPostLast","RHPostLast"];
triggerlabels = cat(3,triggerlabels1,triggerlabels2);

%Locate the EMG recording file corresponding to the specified subject, session, and hand
sessions = strtrim(string(ls(strcat(TMSdir,'\',sid,'\Sessions\')))); %get folder names
sessionnames = sessions(startsWith(sessions,'2'));
sessiontimes = split(sessionnames,'_');
[~,sessionorder] = sort(sessiontimes(:,2)); %sort by time (may be redundant)
sortedsessions = sessionnames(sessionorder); 
filename = char(strcat(TMSdir,'\',sid,"\Sessions\",sortedsessions(isession),"\",sortedsessions(isession),"-emg"));

meps = [];
traces = [];
offset = [];

%Indices of the first and last trials
itrial1 = eval(strcat('trials.',char(triggerlabels(isession,ihand,1)),'(itable)'));
itrial2 = eval(strcat('trials.',char(triggerlabels(isession,ihand,2)),'(itable)'));

%Iterate through trials, loading one at a time using read_eep_trial
for tn = itrial1:itrial2 %get trigger numbers from table
    
    [data,~] = read_eep_trial(filename,tn+1,[-before,after]); %The +1 is very important! Trigger 1 is always "empty"!
    x = data.data;
    fs = data.rate;
    x = x - x(1);
    
    %High-pass Butterworth filter
    [b,a] = butter(1,fc/(fs/2),'high');
    y = filter(b,a,x);
    ltail = ceil(tend*fs); %use "tail" as reference to filter out noise using harmonic analysis
    T = length(x);
    clear s
    s = zeros(length(f),1);

    %Find fundamental frequency within range
    for fi = 1:length(f)
        phase = (0:T-1)'*(1:K)*2*pi*f(fi)/fs;
        S = [sin(phase) cos(phase)]';   % harmonic basis S'
        b = S(:,1:ltail)'\y(1:ltail)';
        s(fi) = std(y'-S'*b);
    end
    [~,fmin]=min(s); %fundamental frequency

    %Remove noise
    phase = (0:T-1)'*(1:K)*2*pi*f(fmin)/fs;
    S = [sin(phase) cos(phase)]';
    b = S(:,1:ltail)'\y(1:ltail)';
    xf = y'-S'*b; %filtered signal
    
    %Detect MEP peaks and valleys
    xsc = (data.time-data.triggers.time*1000)/1000;
    imepp = find((xsc>tmepmin) & (xsc<tmepmax)); %look within specified time range
    [yp,xp] = findpeaks(-xf(imepp),xsc(imepp),"NPeaks",2,"MaxPeakWidth",0.012,"MinPeakProminence",25);  
    if ~isempty(xp) %if peak found
        [~,ip] = max(yp);
        xp = xp(ip);
        yp = yp(ip);
        imepv = find((xsc>=xp) & (xsc<tmepmax));
        [yv,xv] = findpeaks(xf(imepv),xsc(imepv),"NPeaks",2,"MinPeakWidth",0.001,"MinPeakProminence",5);
        if ~isempty(xv) %if peak and valley found
            [~,iv] = max(yv); %pick the lowest valley
            xv = xv(iv);
            yv = yv(iv);
            a = yp+yv; %update amplitude
            meps = [meps;a];
            traces = [traces,xf];
            offset = [offset,data.triggers.time-data.time(1)/1000]; %time of first sample relative to trigger in sec
            % imagesc((1:size(traces,1))/fs-0.5,1:size(traces,2),traces'); hold on %just for checking
        end
    end
end

%Due to a software issue with the TMS/EMG recording system, the recording
%window (relative to the TMS pulse) shifted across trials. This section 
%shifts the samples so that all trials are aligned and overlapping.
t = repmat((1:size(traces,1))/fs,size(traces,2),1)'-repmat(offset,size(traces,1),1); %time
[~,ic] = min(abs(t(:,1))); %index of zero in first trial
tc = t-t(ic,1); %time centered around zero of first trial (the trigger drifts over time)
[~,imins] = min(abs(tc-tc(ic,1)),[],1); %find indices of closest values in subsequent trials to zero in first trial
shift = ic-imins; %number of samples to shift each trial by
alignedtraces = zeros(size(traces));
ttrace = zeros(size(t)); %time for plotting aligned traces
for itrial = 1:length(offset); alignedtraces(:,itrial) = circshift(traces(:,itrial),shift(itrial)); ttrace(:,itrial) = circshift(tc(:,itrial),shift(itrial)); end
alignedtraces([1:max(shift),end-max(shift)+1:end],:)=[]; %remove non-overlapping samples
ttrace([1:max(shift),end-max(shift)+1:end],:)=[]; %remove non-overlapping samples
ttrace = mean(ttrace,2);
traces = alignedtraces;