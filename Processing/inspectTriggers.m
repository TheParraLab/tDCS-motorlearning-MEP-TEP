function EEG = inspectTriggers(session,rawdir,TMSdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Given a EEG recording session, plots TMS pulses detected by TESA and 
%   overlays true TMS triggers from the corresponding EMG recording
%
%   Inputs:
%   session -   recording session identifier, in the same format as the EEG 
%               filename, with subjectid, hand (L/R), and pre/post 
%               stimulation, e.g. 'A33LPre'
%   rawdir -    path to directory with raw EEG recordings
%   TMSdir -    path to directory with EMG recordings. Make sure that the
%               table containing the TMS trial information
%               ('TMStrials.csv') is here.
%   
%   Requires TESA and BVA-iO plugins for EEGLAB. Make sure these are on the
%   MATLAB path!
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,isn] = regexp(session,'\d');
[~,iside] = regexp(session,'\d+.'); %index for hand side (L/R)
itp = regexp(session,'\D');
itp = itp(itp>iside);
tp = session(itp);
hand = session(iside(1));
sn = str2double(session(isn));
sid = strcat('A',num2str(sn));
if strcmp(session(iside(1)),'L') %pick reference channel based on hand
    refchan = 'C4';
elseif strcmp(session(iside(1)),'R')
    refchan = 'C3';
else
    refchan = [];
end
EEG = pop_loadbv(rawdir,strcat(session,'.vhdr')); %load EEG recording
EEG = tesa_findpulse(EEG,refchan,'refract',25,'plots','on');
events = struct2cell(EEG.event);
iEEG = cellfun(@(x) any(strcmp(x,'TMS')),events);
iEEG = find(vertcat(iEEG(7,:,:)));
%Label the detected TMS pulses
for i = 1:length(iEEG)
    text(EEG.event(iEEG(i)).latency,500,num2str(i),'HorizontalAlignment','center')
end

%Load and overlay true TMS trigger events from EMG recordings to check alignment
sessions = strtrim(string(ls(strcat(TMSdir,'\',sid,'\Sessions\')))); %get folder names
sessionnames = sessions(startsWith(sessions,'2'));
sessiontimes = split(sessionnames,'_'); %The TMS recording folders are named after the time of recording
[~,sessionorder] = sort(sessiontimes(:,2)); %sort by time (might be redundant)
sortedsessions = sessionnames(sessionorder);
stp = ["Pre","Post"];
shand = ["L","R"];
isession = find(strcmp(stp,tp));
ihand = find(strcmp(shand,hand));

%Load the triggers from the corresponding EMG recording
filename = char(strcat(TMSdir,"\",sid,"\Sessions\",sortedsessions(isession),"\",sortedsessions(isession),"-emg.trg"));
trg = read_eep_trg(filename);

%Use the TMS trials table to index the corresponding triggers
trialsfile = 'TMStrials.csv'; %file containg TMS trial numbers
trials = readtable(strcat(TMSdir,'\',trialsfile));
triggerlabels1 = ["LHPreFirst","RHPreFirst";"LHPostFirst","RHPostFirst"];
triggerlabels2 = ["LHPreLast","RHPreLast";"LHPostLast","RHPostLast"];
triggerlabels = cat(3,triggerlabels1,triggerlabels2);
itable = find(strcmp(trials.Subject,sid)); %row number matching current subject
itrial1 = eval(strcat('trials.',char(triggerlabels(isession,ihand,1)),'(itable)'));
itrial2 = eval(strcat('trials.',char(triggerlabels(isession,ihand,2)),'(itable)'));

%Scale the times to match EEG recording
tEMG = vertcat(trg(itrial1+1:itrial2+1).time_s)*EEG.srate;
tEMGsc = cumsum(diff(tEMG))+EEG.event(iEEG(1)).latency;
tEMGsc = [EEG.event(iEEG(1)).latency;tEMGsc];
%Label the TMS triggers
plot(tEMGsc,-500,'m*');
for i = 1:length(tEMGsc)
    text(tEMGsc(i),-600,num2str(i),'HorizontalAlignment','center')
end