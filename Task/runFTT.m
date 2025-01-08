%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Starts finger tapping task by calling the SequenceTask2 app
%   Prompts user for subject ID, session number, and sequence number
%   Session 1 = concurrent with tDCS; Session 2 = 1 hour post-tDCS
%   Sequences 1/2/3 = 41324, 23142, 34213
%
%   Copyright (c) 2024 Gavin Hsu. All rights reserved.
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = 'Subject ID: ';
sid = input(prompt,'s');
prompt1 = 'Session number (1/2): ';
session = input(prompt1,'s');
prompt2 = 'Sequence number (1/2/3): ';
sn = input(prompt2,'s');                %Sequence number
while ~strcmp(sn,'1') && ~strcmp(sn,'2') && ~strcmp(sn,'3')
    prompt2 = 'Sequence number (1/2/3): ';
    sn = input(prompt2,'s');
end
filename = strcat('Subject',sid,'_',session,'_S',sn,'.mat');
addpath('BehavioralData')
addpath('Resources')

while isfile(strcat('BehavioralData\',filename)) || isfile(filename)  %Prevent file overwrite
    disp('File already exists. Please try again.')
    prompt = 'Subject ID: ';
    sid = input(prompt,'s');
    prompt1 = 'Session number (1/2): ';
    session = input(prompt1,'s');           %Session number
    prompt2 = 'Sequence number (1/2/3): ';
    sn = input(prompt2,'s');                %Sequence number
    while ~strcmp(sn,'1') && ~strcmp(sn,'2') && ~strcmp(sn,'3')
        prompt2 = 'Sequence number (1/2/3): ';
        sn = input(prompt2,'s');
    end
    filename = strcat('Subject',sid,'_',session,'_S',sn,'.mat');
end

SequenceTask2_App(filename,sn,session)