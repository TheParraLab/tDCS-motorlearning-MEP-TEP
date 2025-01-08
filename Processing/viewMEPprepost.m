%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script loops through all subjects, using plotMEPprepost to
%   visualize MEPs. Make sure that the MEP data are saved in MEPdata.mat
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

load("MEPdata.mat","data")
snmax = 129;
sn = (1:snmax)'; %manually set subject numbers
exclude = []; %subjects to exclude
sn(exclude) = [];

for isub = 1:length(sn)
    plotMEPprepost(sn(isub),data);
    pause
    close all
end