function addeeglab(eeglabpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Adds EEGLAB directories needed for basic functions.
%
%   Input:
%   eeglabpath - path to your EEGLAB directory
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(eeglabpath)
addpath(strcat(eeglabpath,'/functions/sigprocfunc'))
addpath(strcat(eeglabpath,'/functions/popfunc'))
addpath(strcat(eeglabpath,'/functions/guifunc'))
addpath(strcat(eeglabpath,'/functions/adminfunc'))