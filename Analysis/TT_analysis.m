%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script sorts files in folder 'foldername' with file names 
%   containing the prefix 'filetag' and delimeter 'delim' into 'NG' 
%   (3 or more) groups.
%   Group assignments are stored in 'groupfile'
%
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

addpath Resources/
foldername = 'TTData';
filetag = 'TT_data';
delim = '_';
scorevar = 'speed';
groupfile = 'tdcs_groups.csv';
NG = 3;
load('groupnames.mat','groupnames')

%Get file names inside folder containing tag
listing = dir(foldername);
filenames = string(size(listing,1))';
for f = 1:size(listing,1)
    filenames(f) = listing(f).name;
end
filenames = filenames(find(contains(filenames,filetag)));

%Extract subject IDs from file names
subjects = split(filenames,delim);
if length(filenames) > 1
    subjects = reshape(subjects,size(filenames,2),size(subjects,3));
    subjects = split(subjects(:,end)','.');
    subjectids = subjects(:,:,1);
else
    subjects = split(subjects(end),'.');
    subjectids = subjects(1);
end
subjectorder = randperm(length(subjectids)); %shuffle subject order
subjectids = subjectids(subjectorder);

%Load existing group assignments from groupfile
optsTT = detectImportOptions(groupfile);
optsTT.SelectedVariableNames = {'subjectID',scorevar,'group','gn'};
TTT = readtable(groupfile,optsTT);
groupedids = string(table2array(TTT(:,1)));
if isempty(groupedids)
    TTT = convertvars(TTT,{scorevar,'group','gn'},'cell2mat');
end

newids = setdiff(subjectids,groupedids,'stable'); %list of new IDs

%Sort new subjects into groups
nnew = length(newids);
wordspm = zeros(nnew,1);
for isubject = 1:nnew
    filename = strcat(foldername,'/',filetag,delim,newids(isubject),'.csv');
    wordspm(isubject) = get_typingspeed(filename);
    igroup = addtogroup(wordspm(isubject),NG,TTT);
    Cnew = {newids(isubject),wordspm(isubject),groupnames(igroup),igroup};
    disp(strcat("Subject ",num2str(newids(isubject)),": ",num2str(wordspm(isubject))," words per minute - Group ",groupnames(igroup)))
    Tnew = cell2table(Cnew,'VariableNames',{'subjectID',scorevar,'group','gn'});
    TTT = [TTT;Tnew];
end

%Store new assignments and optionally check group balance
writetable(TTT,groupfile)