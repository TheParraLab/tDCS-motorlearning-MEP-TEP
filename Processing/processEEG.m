%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script uses EEGLAB and TESA to processes epoched EEG data through 
%   the following steps to remove TMS and other artifacts:
%   1. Remove TMS artifact
%   2. Remove channel baseline means
%   3. FastICA
%   4. Automatic component classification and rejection using TESA
%   Due to the large total number of trials, we skip manual epoch rejection
%   before ICA.
%
%   This code was originally created by Davide Bonfanti and has been
%   modified by Gavin Hsu.
%   It is based on suggestions by the creators of TESA. For more info:
%   https://www.sciencedirect.com/science/article/pii/S1053811916305845 
%   https://nigelrogasch.gitbook.io/tesa-user-manual/
%
%   The TESA plugin for EEGLAB is required to run this script (tested with
%   1.1.1). Make sure it's on the MATLAB path!
%   FastICA is required to run this script (tested with 2.5). Make sure 
%   it's on the MATLAB path!
%   For more info: https://research.ics.aalto.fi/ica/fastica/
% 
%   Copyright (c) 2024 Davide Bonfanti
%   Modified (c) 2024 Gavin Hsu, Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

maindir = input('Enter the path to your main directory: ',"s");
addpath(maindir)
eeglabpath = input('Enter the path to your EEGLAB directory: ',"s"); %Path to EEGLAB
addeeglab(eeglabpath)
edir = input('Enter the path to your epoched data: ',"s");
pdir = strcat(edir,'\Processed\');  %This is where the processed data will be stored
if ~exist(pdir,'dir')
    mkdir(pdir)
end

files = strtrim(string(ls(strcat(edir,'*.set')))); %list of epoched data files
sns = zeros(length(files),1); %list of subject ID numbers
%Get list of unique subject ID numbers
for i = 1:length(files)
    fileparts = split(files(i),'.');
    session = char(fileparts(1)); %extract session name
    [~,isn] = regexp(session,'\d'); %extract indices for subject number
    sns(i) = str2double(strcat(session(isn)));
end
usn = unique(sns); %list of subjects found

%Combine pre/post tDCS epochs so that they are processed consistently
for isn = 1:length(usn)
    combineEEG(usn(isn),edir,pdir)
end

cd(pdir)

%%  1.Remove TMS artifact 
files = strtrim(string(ls(strcat(pdir,'*_c.set')))); %list of combined epoch datasets

for k=1:length(files)
    %Check for existing processed datasets
    if exist([files{k}(1:end-4),'_RTA.set'],'file')
        disp(strcat("Skipping ",files{k}(1:end-4),'_RTA.set'))
    else
        ALLEEG = [];
        EEG = pop_loadset('filename',files{k},'filepath',pdir);
        EEG = eeg_checkset(EEG);
        [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);

        %Remove samples around TMS pulse (time = 0)
        EEG = pop_tesa_removedata(EEG,[-3,10],[-500,-150]); %remove -3 to 10 ms, referencing -500 to -150ms to replace the samples
        EEG = eeg_checkset(EEG);
        pop_newset(ALLEEG, EEG, 1,'savenew',[files{k}(1:end-4),'_RTA.set'],'gui','off');
        clear EEG
    end
end


%% 2. Remove channel baseline means
files = strtrim(string(ls(strcat(pdir,'*RTA.set')))); %list of processed datasets from the previous step

for k=1:length(files)
    %Check for existing processed datasets
    if exist([files{k}(1:end-4),'_DM.set'],'file')
        disp(strcat("Skipping ",files{k}(1:end-4),'_DM.set'))
    else
        ALLEEG = [];
        EEG = pop_loadset('filename',files{k},'filepath',pdir);
        EEG = eeg_checkset(EEG);
        [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);

        %Demean
        EEG = pop_rmbase(EEG,[-500,-400]); %take the baseline from -500 to -400ms
        EEG = eeg_checkset(EEG);
        pop_newset(ALLEEG,EEG,1,'savenew',[files{k}(1:end-4),'_DM.set'],'gui','off');
        clear EEG
    end
end

%% 3. FastICA
%Due to the large number of subjects and trials, we do not manually inspect
%every individual epoch prior to ICA.
files = strtrim(string(ls(strcat(pdir,'*_DM.set')))); %list of processed datasets from the previous step

for k=1:length(files)
    %Check for existing processed datasets
    if exist([files{k}(1:end-4),'_fastICA1.set'],'file')
        disp(strcat("Skipping ",files{k}(1:end-4),'_fastICA1.set'))
    else
        EEG = pop_loadset('filename',files{k},'filepath',pdir);
        EEG = eeg_checkset(EEG);
        [ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG,0);
        
        %Run FastICA
        EEG = pop_tesa_fastica(EEG,'approach','symm','g','tanh','stabilization','on'); %symmetric approach, Gauss contrast function, with stabilized FastICA
        EEG = eeg_checkset(EEG);
        pop_newset(ALLEEG,EEG,1,'savenew',[files{k}(1:end-4),'_fastICA1.set'],'gui','off');
    end
end

%% 4. Automatic component classification and rejection using TESA
%Again, due to the large amount of subjects and trials, we do not manually
%reject components and instead rely on TESA.
files = strtrim(string(ls(strcat(pdir,'*_fastICA1.set')))); %list of processed datasets from the previous step

for k=1:length(files)
    %Check for existing processed datasets
    if exist([files{k}(1:end-4),'_TESA.set'],'file')
        disp(strcat("Skipping ",files{k}(1:end-4),'_TESA.set'))
    else
        EEG = pop_loadset('filename',files{k},'filepath',pdir);
        EEG = eeg_checkset(EEG);

        %TESA classification
        EEG = tesa_compselect(EEG,'compCheck','off','remove','off','plotTimeX',[-200,499]);
        artcomps = find(EEG.icaCompClass.TESA1.compClass > 1); %Components identified as artifacts from muscle, blinks, eye movements, electrode noise, sensory, etc.
        %Automatically identify step artifacts with slope above a threshold
        stepcomps = findstepIC(EEG);
        artcomps = [artcomps,stepcomps'];
    
        EEG = pop_subcomp(EEG,artcomps,0,0);
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG,'filepath',pdir,'filename',[files{k}(1:end-4),'_TESA.set']);
        clear EEG
    end
end

%% Optional: plot results
% filesfx = '*_TESA.set';
% files = strtrim(string(ls(strcat(pdir,filesfx)))); %list of epoched data files
% sns = zeros(length(files),1); %list of subject numbers
% for i = 1:length(files)
%     fileparts = split(files(i),'.');
%     session = char(fileparts(1)); %extract session name
%     [isn1,isn2] = regexp(session,'\d+\D'); %extract indices for subject number
%     sns(i) = str2double(strcat(session(isn1(1):isn2(1)-1)));
% end
% usn = unique(sns);
% for isn = 1:length(usn)
%     plotTEP(usn(isn),filesfx)
%     movegui('center')
%     pause
% end
