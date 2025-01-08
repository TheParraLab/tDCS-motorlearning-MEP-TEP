%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script loads processed TEP data and stores them into a .mat file
%   for analysis. It prompts the user to enter the directories for EEGLAB
%   and the processed data.
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

eeglabpath = input('Enter the path to your EEGLAB directory: ',"s"); %Path to EEGLAB
addeeglab(eeglabpath)
pdir = input('Enter the path to your processed data: ',"s");
filesfx = '*_TESA.set'; %file suffix of processed data
files = strtrim(string(ls(strcat(pdir,filesfx)))); %list of processed data files
sns = zeros(length(files),1); %list of subject numbers

%Get unique subject numbers from the file list
for i = 1:length(files)
    fileparts = split(files(i),'.');
    session = char(fileparts(1)); %extract session name
    [isn1,isn2] = regexp(session,'\d+\D'); %extract indices for subject number
    sns(i) = str2double(strcat(session(isn1(1):isn2(1)-1)));
end
usn = unique(sns);

hands = ["L","R"];
sessions = ["Pre","Post"];
TEPdata = struct; %struct where all TEP data will be stored

%Load data and store into struct
for isn = 1:length(usn)
    sn = usn(isn);
    TEPdata(isn).sn = sn;
    sid = strcat('A',strtrim(string(num2str(sn(:)))));
    filenames = strtrim(string(ls(strcat(pdir,sid,filesfx))));
    for ihand = 1:length(hands)
        EEG = pop_loadset('filename',char(filenames(ihand)),'filepath',pdir);
        istart = [1;EEG.ntrials1+1]; %index for first trial in each session
        iend = [EEG.ntrials1;EEG.ntrials2+EEG.ntrials1]; %index for last trial in each session
        for isession = 1:length(sessions)
            datamt = median(EEG.data(:,:,istart(isession):iend(isession)),3);
            eval(strcat('TEPdata(isn).datamt',hands(ihand),sessions(isession),' = datamt;'))
        end
    end
end
%Save struct into .mat file, along with other info
m = matfile('TEPdata.mat','Writable',true);
m.TEPdata = TEPdata;
m.t = EEG.times; %time stamps of EEG samples
m.fs = EEG.srate; %sampling rate
m.nsamples = size(EEG.times,2); %number of samples in each epoch
m.chanlocs = EEG.chanlocs; %EEG channels
m.nchan = EEG.nbchan; %number of EEG channels
disp('TEP data saved')
clear m