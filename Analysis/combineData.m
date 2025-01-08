%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script is used to combine behavioral, MEP, and TEP data into a
%   structure array, stored in 'combinedData.mat'
%   The .mat file is loaded in figsandstats.m for quantitative analysis and
%   generating figures
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

maindir = input('Enter the path to your main directory: ',"s"); %Path to main directory
Bfolder = 'BehavioralData/'; %Folder containing behavioral data (under main directory)
TTfolder = 'TTData'; %Folder containing typing test data (under main directory)

cd(maindir)
addpath Resources/
filetag = 'TT_data'; %Filename tag for typing test data
delim = '_';

%Get file names inside folder that contain the tag
listing = dir(TTfolder);
filenames = string(size(listing,1))';
for f = 1:size(listing,1)
    filenames(f) = listing(f).name;
end
ifiles = find(contains(filenames,filetag));
filenames = filenames(ifiles);

ns = length(filenames); %Total number of subjects

% %For testing with randomized groups: extract subject IDs from file names
% subjects = split(filenames,delim);
% if length(filenames) > 1
%     subjects = reshape(subjects,size(filenames,2),size(subjects,3));
%     subjects = split(subjects(:,end)','.');
%     subjectids = subjects(:,:,1);
% else
%     subjects = split(subjects(end),'.');
%     subjectids = subjects(1);
% end
% sids = str2double(subjectids(:));
% [sids,isids] = sort(sids); %Subject IDs in ascending order
% % groupsids = sids(randperm(ns)); %Just for testing!

%Load group info
opts = detectImportOptions("tdcs_groups.csv");
opts.SelectedVariableNames = ["subjectID","gn"];
GT = readtable("tdcs_groups.csv",opts); %Grouping table
sns = GT.subjectID; %Subject IDs corresponding to behavioral data
grouporder = [3,1,2]; %group numbers corresponding to 0, 4, 6 mA, respectively

%Load demographics info
SAS = readtable('demographics.csv');

ng = 3; %Number of groups
S = struct;
nsg = floor(ns/ng); %Number of subjects per group
ntrials = 36; %Number of trials per task
nstrials = 10; %Number of trials to include in "early" and "saturated"
nsec = 4; %Number of sections/tasks
siit = 1:16;    %Data for first 16 subjects contain an invalid trial
sessions = ["Pre","Post"];

%Load MEP data
MEPdata = load('MEPdata.mat','data');
MEPdata = MEPdata.data;
MEPsn = vertcat(MEPdata.sn); %Subject IDs corresponding to MEP data
% MEPpostpre = vertcat(MEPdata.postpre); %MEP post/pre ratios: left hand in column 1 and right hand in column
MEPcell = struct2cell(MEPdata);
allfn = fieldnames(MEPdata);

%Load TEP data
TEPfile = load('TEPdata.mat','TEPdata','chanlocs','t','nchan','nsamples');
EEGns = TEPfile.nsamples; %Number of samples per EEG epoch
nchan = TEPfile.nchan; %Number of channels in EEG
EEGt = TEPfile.t; %Time for EEG recording
chanlocs = TEPfile.chanlocs; %EEG channel locations
TEPdata = TEPfile.TEPdata;
TEPsn = vertcat(TEPdata.sn); %Subject IDs corresponding to TEP data
pt = [33,45,55,100,175]'; %Peak times to evaluate
hw = 5; %Half width of time window around peak times

for ig = 1:ng
    % S(ig).sid = groupsids((ig-1)*nsg+1:(ig-1)*nsg+nsg)'; %Just for testing!
    S(ig).sid = GT.subjectID(GT.gn==grouporder(ig)); %Get subject IDs belonging to the current group
    S(ig).n = nsg; %Number of subjects per group
    %Trialwise measurements: each column represents one subject and each row represents one trial
    S(ig).ts_all = zeros(ntrials,nsg,nsec);   %Tapping speed for all trials, subjects, sections
    S(ig).a_all = zeros(ntrials,nsg,nsec);    %Accuracy for all trials, subjects, sections
    S(ig).ncs_all = zeros(ntrials,nsg,nsec);   %Number of correct sequences for all trials, subjects, sections
    S(ig).r = zeros(3,nsg);    %Sensation rating
    S(ig).nr = zeros(nsec,nsg);   %Number of trials with no response in each section
    S(ig).lp = zeros(nsec,nsg);   %Number of trials with long pauses in each section
    S(ig).wordspm = zeros(nsg,1);   %Typing speed
    S(ig).iMEP = NaN(nsg,1);
    S(ig).MPPL = NaN(nsg,1);    %MEP post/pre ratio for left hand
    S(ig).MPPR = NaN(nsg,1);    %MEP post/pre ratio for right hand
    S(ig).TEPL = NaN(nsg,nchan,EEGns,2);    %TEP for left hand - subject, channel, samples, pre/post
    S(ig).TEPR = NaN(nsg,nchan,EEGns,2);    %MEP for right hand - subject, channel, samples, pre/post
    S(ig).nkp = zeros(nsec,nsg);    %Number of key presses throughout task
    S(ig).TOD = zeros(nsg,1);   %Time of day (0 = AM, 1 = PM)
    for ihand = 1:length(hands)
        for isession = 1:length(sessions)
            eval(strcat('S(ig).MEP',hands(ihand),sessions(isession),'=NaN(nsg,1);'))
        end
    end

    %Calculate performance
    for isubject = 1:nsg %isubject is the subject index within the group, NOT the subject ID!
        disp(strcat("Gathering data for ",num2str(S(ig).sid(isubject)),"..."))
        filename = strcat(TTfolder,'/',filetag,delim,num2str(S(ig).sid(isubject)),'.csv');
        S(ig).wordspm(isubject) = get_typingspeed(filename);
        [S(ig).ts_all(:,isubject,:),S(ig).a_all(:,isubject,:),S(ig).ncs_all(:,isubject,:),S(ig).r(:,isubject),~,~,S(ig).nr(:,isubject),S(ig).lp(:,isubject),S(ig).nkp(:,isubject)] = get_FTT(S(ig).sid(isubject),Bfolder);
        iMEP = find(MEPsn==S(ig).sid(isubject),1);
        if ~isempty(iMEP)
            S(ig).iMEP(isubject) = iMEP;
            for ihand = 1:length(hands)
                for isession = 1:length(sessions)
                    fn = strcat('meps',hands(ihand),'H',sessions(isession)); %Field name to search for
                    ifn = find(strcmp(allfn,fn)); %Field index
                    MEPfield = MEPcell(ifn,:,:);
                    if size(MEPfield{:,:,iMEP},1)>30
                        eval(strcat('S(ig).MEP',hands(ihand),sessions(isession),'(isubject)=squeeze(cellfun(@median,MEPfield(:,:,iMEP)));'))
                    end
                end
            end
            S(ig).MPPL(isubject) = S(ig).MEPLPost(isubject)/S(ig).MEPLPre(isubject);
            S(ig).MPPR(isubject) = S(ig).MEPRPost(isubject)/S(ig).MEPRPre(isubject);
        end
        iTEP = find(TEPsn==S(ig).sid(isubject),1);
        if ~isempty(iTEP)
            for ihand = 1:length(hands)
                for isession = 1:length(sessions)
                    eval(strcat('S(ig).TEP',hands(ihand),'(isubject,:,:,',num2str(isession),')=TEPdata(iTEP).datamt',hands(ihand),sessions(isession),';'))
                end
            end
        end
        filetime = datetime(dir(filename).datenum,'ConvertFrom','datenum','format','H');
        S(ig).TOD(isubject) = double(hms(filetime)>=12);
        idemo = find(table2array(SAS(:,1))==S(ig).sid(isubject));
        S(ig).age(isubject) = table2array(SAS(idemo,2));
        S(ig).sex(isubject) = table2array(SAS(idemo,3));
    end
    S(ig).tTEP = EEGt;
    
    %Omit the invalid trials
    iit = zeros(size(S(ig).ncs_all)); %logical indices of invalid trials
    [~,isiit,~] = intersect(S(ig).sid,siit); %within-group subject index that match subject IDs with invalid trials
    iit(1,isiit,:) = 1; %mark the invalid trials
    iit = logical(iit);
    S(ig).ts_all(iit) = NaN;
    S(ig).a_all(iit) = NaN;
    S(ig).ncs_all(iit) = NaN;
    niit = ~iit;
    [~,ifvt] = max(niit(:,:,1),[],1); %index of the first valid trials
    
    %N Correct Sequences: Mean "Saturated" (in the Last nstrials)
    S(ig).ncs_sat = S(ig).ncs_all(ntrials-(nstrials-1):ntrials,:,:);
    S(ig).ncs_sat_m = mean(S(ig).ncs_sat,1,'omitnan');
    
    %N Correct Sequences: Mean "Early" (in the First nstrials)
    S(ig).ncs_ear = S(ig).ncs_all(ifvt:ifvt+nstrials,:,:);
    S(ig).ncs_ear_m = mean(S(ig).ncs_ear,1,'omitnan');
    
    %N Correct Sequences: Mean Across All Trials
    S(ig).mncst = mean(S(ig).ncs_all,1,'omitnan');
    S(ig).mncst = reshape(S(ig).mncst,S(ig).n,4);

    %N Correct Sequences: Mean Across Subjects
    S(ig).mncss = mean(S(ig).ncs_all,2,'omitnan');
    nvalid = S(ig).n-sum(isnan(S(ig).ncs_all),2); %Number of valid samples
    S(ig).encss = std(S(ig).ncs_all,0,2,'omitnan')./sqrt(nvalid); %SEM
    
    %Accuracy: Mean Across All Trials
    S(ig).a_m = mean(S(ig).a_all,1,'omitnan');
    
    %Change in NCS
    S(ig).ncs_start = zeros(1,nsg);
    for isubject = 1:nsg
        temp = S(ig).ncs_all(:,isubject,:);
        temp = temp(~isnan(temp(:,:,1)));
        S(ig).ncs_start(isubject) = temp(find(temp,1));
    end
    S(ig).ncs_diff = (S(ig).ncs_sat_m - S(ig).ncs_start);
    S(ig).ncs_diff = reshape(S(ig).ncs_diff,nsg,nsec);
    
    %Sensation
    S(ig).mr = mean(S(ig).r,2);
    S(ig).er = std(S(ig).r,0,2)/sqrt(size(S(ig).r,2)); %SEM

    %Tapping Speed Mean Across All Trials
    S(ig).mtst = mean(S(ig).ts_all,1,'omitnan');
    S(ig).mtst = reshape(S(ig).mtst,S(ig).n,4);

    %Check for normal distribution (returns 1 if not normal)
    S(ig).normal = zeros(4,1);
    for itask = 1:4
        S(ig).normal(itask) = lillietest(S(ig).mncst(:,itask));
    end
    S(ig).normalMPPL = lillietest(S(ig).MPPL);
    S(ig).normalMPPR = lillietest(S(ig).MPPR);

    %TEP values at specified peak times - subject, channel, left/right, pre/post
    for ip = 1:length(pt)
        for ihand = 1:length(hands)
            eval(strcat('S(ig).P',num2str(pt(ip)),'(:,:,',num2str(ihand),',:)=squeeze(mean(S(ig).TEP',hands(ihand),'(:,:,find(EEGt<',num2str(pt(ip)),'+hw&EEGt>',num2str(pt(ip)),'-hw),:),3));'))
        end
    end
end

clear MEPdata
clear TEPdata
clear TEPfile

disp('Saving .mat file...')
m = matfile('combinedData.mat','Writable',true);
m.S = S;
disp('Done')