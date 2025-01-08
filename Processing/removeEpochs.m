function removeEpochs(session,ireject,edir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Given a EEG recording session, removes bad trials from epoched EEG data
%   Warning: this overwrites the original epoched dataset!
%
%   Inputs:
%   session -   recording session identifier, in the same format as the EEG 
%               filename, with subjectid, hand (L/R), and pre/post 
%               stimulation, e.g. 'A33LPre'
%   ireject -   trial indices of trials to remove, e.g. ireject = [61:82];
%   edir -      path to directory with epoched EEG data
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG = pop_loadset([edir,session,'.set']);
itrials = zeros(EEG.trials,1);
itrials(ireject) = 1;
EEG2 = pop_rejepoch(EEG,itrials,0); 
pop_saveset(EEG2,'filepath',edir,'filename',session);