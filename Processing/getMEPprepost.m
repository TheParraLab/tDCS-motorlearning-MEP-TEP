%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script processes MEP data for each subject, storing MEP amplitudes
%   and aligned EMG traces.
%   Before running, make sure that trigger numbers were entered into the
%   table 'TMStrials.csv'
%   Pre/post and L/R data for each subject are stored in a struct 'newdata'
%   If not all subjects are processed at once, any newly processed data
%   must be concatenated with previous samples and saved manually, e.g.:
%   data = [data,newdata];
%   m = matfile('MEPdata.mat','Writable',true);
%   m.data = data;
%
%   Inputs:
%   TMSdir -    directory where raw EMG recordings are stored
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

TMSdir = input('Enter the path to your raw EMG data: ',"s");

%Set the subjects to process
snmax = 129;
sn = (1:snmax)';
exclude = [1:5,11,12,20,42,46]; %subjects to exclude
sn(exclude) = [];

sessions = ["Pre","Post"];
hands = ["LH","RH"];
newdata = struct;

%Load existing data to avoid overwriting
if exist('MEPdata.mat','file')
    load("MEPdata.mat","data")
    sns = vertcat(data.sn);
else
    sns = [];
end
isubnew = 1;

%Go through all subjects and store MEP data into newdata
isnormal = NaN(length(sn),length(hands),length(sessions)); %Lillieforst test decision: subject, session, hand
for isub = 1:length(sn)
    %Check for existing data
    isn = find(sns==sn(isub),1); 
    if ~isempty(isn)
        disp(strcat("Subject ",num2str(sn(isub))," already recorded!"))
    else
        for isession = 1:2
            meps = [];
            traces = [];
            ttrace = [];
            for ihand = 1:2
                %Get MEPs from EMG recording
                [meps,traces,ttrace,fs] = getMEPs(sn(isub),isession,ihand,TMSdir);
                %Find outliers
                logpow = log(meps.^2);
                qr3 = 3*iqr(meps);
                medmep = median(meps);
                outliers = meps < (medmep-qr3) | meps > (medmep+qr3); %3 quartiles from median
                disp(strcat(num2str(sum(outliers))," outliers found"))
                %Store into struct
                eval(strcat('newdata(isubnew).meps',hands(ihand),sessions(isession),' = meps(~outliers);'))
                eval(strcat('newdata(isubnew).traces',hands(ihand),sessions(isession),' = traces;'))
                eval(strcat('newdata(isubnew).t',hands(ihand),sessions(isession),' = ttrace;'))  
                eval(strcat('newdata(isubnew).outliers',hands(ihand),sessions(isession),'= find(outliers);'))
                eval(strcat('newdata(isubnew).noutliers',hands(ihand),sessions(isession),' = sum(outliers);'))
                %Check for normal distribution
                if length(meps(~outliers)) >= 4
                    isnormal(isubnew,ihand,isession) = lillietest(meps(~outliers));
                end
            end
        end
        newdata(isubnew).fs = fs;
        newdata(isubnew).sn = sn(isub);
        %Calculate post/pre MEP ratio
        for ihand = 1:2
            eval(strcat('newdata(isubnew).postpre(ihand) = median(newdata(isubnew).meps',hands(ihand),sessions(2),')/median(newdata(isubnew).meps',hands(ihand),sessions(1),');'))
        end
        disp(strcat("Subject ",num2str(sn(isub))," done"))
        isubnew = isubnew + 1;
    end
end
disp("Remember to save the new data!")

nnnormal = sum(isnormal,[2,3]); %Number of non-normal distributions for each subject