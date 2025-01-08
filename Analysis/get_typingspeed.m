function wordspm = get_typingspeed(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculates typing speed (words per minute) for a given typing test
%   results .csv file. The typing test files do not contain keystroke 
%   information, since timing is only recorded when the correct key is 
%   pressed.
%
%   Input:
%   filename -  name of the .csv file containing the typing test timing
%
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numwords = 62; %Total number of words in the passage. 

optsTT = detectImportOptions(filename);
optsTT.SelectedVariableNames = {'RT','ntask'};
TTT = readmatrix(filename,optsTT);
TTT = TTT(find(TTT(:,2)==1),:);
typingtime = sum(TTT(:,1));
wordspm = numwords/(typingtime/60);