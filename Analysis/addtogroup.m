function igroup = addtogroup(score,NG,groupinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sorts a subject into a group 'igroup' based on their 'score', given 
%   prior group assignments. The algorithm minimizes variance across groups
%
%   Inputs:
%   score -     The parameter used to sort subjects by
%   NG -        Number of groups
%   groupinfo - Table containing previous group assignments, with
%               previous subject IDs in column 1, corresponding scores in
%               column 2, and group assignments and column 3.
%
%   Output:
%   igroup -    Index of the group for the new subject
%
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(groupinfo,1) < NG %The first sample in each group is sorted sequentially
    igroup = mod(size(groupinfo,1),NG)+1;
else
    scorevar = groupinfo.Properties.VariableNames{2};
    x = zeros(NG,1);
    y = x;
    z = x;
    groupsizes = zeros(NG,1);
    for ig = 1:NG
        irows = find(groupinfo.gn == ig);
        groupsizes(ig) = length(irows);
        eval(strcat('score',num2str(ig),' = groupinfo.',scorevar,'(irows);'))
        eval(strcat('score',num2str(ig),'temp = [score',num2str(ig),';',num2str(score),'];')) %Temporarily add new sample to each group
        x(ig) = mean(eval(strcat('score',num2str(ig),'temp')));
        y(ig) = mean(eval(strcat('score',num2str(ig))));
    end
    for ig = 1:NG
        iy = [1:ig-1,ig+1:length(y)];
        z(ig) = var([y(iy);x(ig)]);
    end
    [~,grouprank] = sort(z);
    [nmin,lastgroup] = min(groupsizes);
    if mod(size(groupinfo,1),NG) == NG - 1
        igroup = lastgroup;
    else
        imin = find(groupsizes == nmin); %Smallest group(s). Includes all groups is group sizes are equal
        minrank = intersect(grouprank,imin,'stable'); %Pick lowest variance from among smallest group(s)
        igroup = minrank(1);
    end
end