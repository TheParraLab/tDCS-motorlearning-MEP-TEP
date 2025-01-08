function combineEEG(sn,edir,pdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Combines pre/post tDCS EEG epochs
%
%   Inputs:
%   sn -      Subject ID number, e.g. 33
%   edir -    path to directory with raw EEG recordings
%   pdir -    path to directory with epoched EEG data
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(edir)
addpath(pdir)
sid = strcat('A',strtrim(string(num2str(sn))));
hands = ["L","R"];
sessions = ["Pre","Post"];
filesfx = '.set';
files = strtrim(string(ls(strcat(pdir,'*_c.set'))));

if ~any(strcmp(files,char(strcat(sid,hands(1),'_c.set'))))
    for ihand = 1:2
        for isession = 1:2
            filename = char(strcat(sid,hands(ihand),sessions(isession),filesfx));
            eval(strcat('EEG',num2str(isession),'= pop_loadset(''filename'',filename,''filepath'',edir);'))
        end
        EEG = pop_mergeset(EEG1,EEG2);
        EEG.ntrials1 = EEG1.trials;
        EEG.ntrials2 = EEG2.trials;
        pop_saveset(EEG,'filepath',pdir,'filename',char(strcat(sid,hands(ihand),'_c.set')));
    end
else
    disp(strcat(sid," already combined"))
end