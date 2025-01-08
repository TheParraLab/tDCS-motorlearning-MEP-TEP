function checkTMSTriggers(sn,TMSdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Checks for extra triggers in trigger recordings.
%   TMS pulses were delivered with at least 5s intervals in between, so any
%   triggers with intervals any shorter are likely false.
%
%   Inputs:
%   sn -        Subject ID number, e.g. 33
%   TMSdir -    Directory where raw EMG recordings are stored
%
%   Requires Libeep MATLAB importer (tested with 3.3.177)
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sid = strcat('A',num2str(sn));
found = 0;

%Sort folders to locate the files in chronological order
sessions = strtrim(string(ls(strcat(TMSdir,'\',sid,'\Sessions\')))); %get folder names
sessionnames = sessions(startsWith(sessions,'2'));
sessiontimes = split(sessionnames,'_');
[~,sessionorder] = sort(sessiontimes(:,2)); %sort by time (may be redundant)
sortedsessions = sessionnames(sessionorder); 

%Open the trigger files
for isession = 1:2
    filename = char(strcat(TMSdir,'\',sid,"\Sessions\",sortedsessions(isession),"\",sortedsessions(isession),"-emg.trg"));
    [trg] = read_eep_trg(filename);
    if trg(3).offset - trg(2).offset < 10000 %Look for offset difference less than 10000 between triggers 2 and 3
        disp(strcat("Extra triggers found in ",filename))
        found = 1;
    end
end

if found == 0
    disp("No extra triggers found")
end