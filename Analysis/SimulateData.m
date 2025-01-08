%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script simulates an intensity effect by altering a copy of 
%   existing data.
%   Simulated data are used to check whether the analysis code is able to
%   correctly detect an effect or lack thereof.
%
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath Resources/
Sfolder = 'SimData/'; %Folder to deposit simulated data in
TTfolder = 'TTData'; %Typing test data are stored here
filetag = 'TT_data'; %Filename prefix for typing test data
delim = '_';
copyfile('BehavioralData/Subject*.mat',Sfolder)

%Get file names inside folder that contain the tag
listing = dir(TTfolder);
filenames = string(size(listing,1))';
for f = 1:size(listing,1)
    filenames(f) = listing(f).name;
end
ifiles = find(contains(filenames,filetag));
filenames = filenames(ifiles);
ns = length(filenames); %Total number of subjects
ng = 3; %Number of groups
nsg = floor(ns/ng); %Number of subjects per group

%Load group info
opts = detectImportOptions("tdcs_groups.csv");
opts.SelectedVariableNames = ["subjectID","gn"];
GT = readtable("tdcs_groups.csv",opts); %Grouping table
grouporder = [3,1,2]; %group numbers corresponding to 0, 4, 6mA, respectively

sidsog = zeros(ns,1); %Subject IDs ordered by group
for ig = 1:ng
    sidsog((ig-1)*nsg+1:(ig-1)*nsg+nsg) = GT.subjectID(GT.gn==grouporder(ig)); %Get subject IDs belonging to the current group
end
current = [zeros(nsg,1);4*ones(nsg,1);6*ones(nsg,1)]; %Current

ntrials = 36;   %Number of trials per section
nsections = 4;  %Number of sections/tasks

for j = 1:ns
    ises = 1;   %Session counter: two sessions total: 1 and 2
    iseq = 1;   %Sequence counter: only one for session 1, three different sequences for session 2
    for m = 1:nsections
        clear log
        session = num2str(ises);
        sn = num2str(iseq);
        filename = strcat('Subject',num2str(sidsog(j)),'_',session,'_S',sn,'.mat');
        load(strcat(Sfolder,filename),'log')
    
        if ~exist('log','var')
            disp(filename)
        else
            for i = 1:ntrials    
                for ikey = 1:size(log.(strcat('trial',num2str(i))),1)
                    %Insert errors in the keystrokes at different rates depending on current intensity
                    if rand(1)<0.01*current(j)
                        log.(strcat('trial',num2str(i)))(ikey,1) = "x";
                    end
                end
            end
            save(strcat(Sfolder,filename),'log','-append')
        end
        if ises ~= 1
            iseq = iseq + 1;
        end
        if ises == 1
            ises = ises + 1;
        end
    end
end