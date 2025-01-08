function T = saveEpochs(sid,rawdir,edir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Extracts and saves epochs from TMSEEG recordings corresponding to the
%   subject sid from the raw folder rawdir and stores epochs into the 
%   folder edir
%   Requires TESA and BVA-iO plugins for EEGLAB
%   Raw files have filenames containing the subject ID, hand  (L/R), and 
%   session (pre/post). For example: A33LPost
%   There should be 3 files for each recording: .eeg, .vhdr, and .vmrk
%
%   Inputs:
%   sid -       Subject ID, e.g. 33
%   rawdir -    path to directory with raw EEG recordings
%   edir -      path to directory with epoched EEG data
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First, check existing epoched data to prevent overwriting
files = strtrim(string(ls(strcat(edir,'\*.set')))); %list of epoched data files
sns = zeros(length(files),1); %list of subject numbers
%Get subject numbers from the file list
for i = 1:length(files)
    fileparts = split(files(i),'.');
    session = char(fileparts(1)); %extract session name (subj ID + L/R + pre/post)
    [~,isn] = regexp(session,'\d'); %extract indices for subject number
    sns(i) = str2double(strcat(session(isn)));
end
%Skip this subject if data already exist
if ~isempty(find(sns==sid,1))
    disp(strcat("Epoched data files already exist for subject ",num2str(sid),"!"))
    T = [];
    return
end

%Look for files in rawdir with matching sid
clear files
clear sns
files = strtrim(string(ls(strcat(rawdir,'\*.vhdr')))); %list of .vhdr files
sns = zeros(length(files),1); %list of subject numbers
%Get subject numbers from the file list
for i = 1:length(files)
    fileparts = split(files(i),'.');
    session = char(fileparts(1)); %extract session name
    [~,isn] = regexp(session,'\d'); %extract indices for the subject number
    sns(i) = str2double(strcat(session(isn)));
end
iuse = find(sns==sid); %find indices matching sid
sn = zeros(length(iuse),1); %redundant list of subject numbers for QC
nstims = zeros(length(iuse),1); %number of TMS stims detected
sessions = strings(length(iuse),1); %session name: subjectid+hand+timepoint

%Extract and save epochs from the matching files
for ifile = 1:length(iuse)
    file = files(iuse(ifile));
    fileparts = split(file,'.');
    session = char(fileparts(1)); %extract session name
    [~,isn] = regexp(session,'\d');
    [~,iside] = regexp(session,'\d+.'); %index for hand side (L/R)

    %Check which hand was stimulated for the recording and pick reference channel based on hand
    if strcmp(session(iside(1)),'L')
        refchan = 'C4';
    elseif strcmp(session(iside(1)),'R')
        refchan = 'C3';
    else
        refchan = [];
    end

    disp(session)
    %Use TESA to find TMS pulses
    EEG = pop_loadbv(rawdir,file); %load EEG recording
    EEG = tesa_findpulse(EEG,refchan,'refract',25,'plots','off'); %find pulses
    events = struct2cell(EEG.event);
    %Look for 'Buffer Overflow' events and delete them
    iBO = cellfun(@(x) any(strcmp(x,'Buffer Overflow')),events);
    iBO = find(vertcat(iBO(7,:,:)));
    if ~isempty(iBO)
        tBO = EEG.event(iBO).latency;
        irm = find(vertcat(EEG.event.latency)==tBO);
        EEG.event(irm) = [];
        disp('Buffer Overflow found!')
    end
    %Extract and save 1-second long epochs
    EEG = pop_epoch(EEG,'TMS',[-0.5,0.5]);
    pop_saveset(EEG,'filepath',edir,'filename',char(fileparts(1)));

    sn(ifile) = str2double(session(isn));
    sessions(ifile) = session;
    nstims(ifile) = size(EEG.data,3);
end
T = table(sn,sessions,nstims); %return number of stimulations detected for each recording, for inspection