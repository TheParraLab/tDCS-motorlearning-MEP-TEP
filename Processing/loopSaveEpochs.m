%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script loops through multiple subjects to extract and save epochs 
%   from raw EEG recordings. Epoched data are stored into the
%   user-specified folder.
%   Requires TESA and BVA-iO plugins for EEGLAB
%   The table TT shows the number of TMS stimuli found in each recording.
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eeglabpath = input('Enter the path to your EEGLAB directory: ',"s");
addeeglab(eeglabpath)
addpath(strcat(eeglabpath,'\plugins\TESA1.1.1\'))
addpath(strcat(eeglabpath,'\plugins\bva-io1.73\'))
rawdir = input('Enter the path to your raw data: ',"s");
edir = input('Enter the path to your epoched data: ',"s");
snmax = 129; %Upper limit of subject ID numbers to process

TT = [];
for sn = 1:snmax
    T = saveEpochs(sn,rawdir,edir);
    TT = [TT;T]; %This table summarizes the number of stimuli found in each recording
end                         