function plotTEP(sn,filesfx,pdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plots processed TEP data for a specified subject, averaged across 
%   trials and channels.
%
%   Inputs:
%   sn -        Subject ID number, e.g. 33
%   filesfx -   file suffix for the processing step to visualize, e.g.
%               "TESA" or "fastICA1" (see processEEG.m)
%   pdir -      path to directory with processed EEG data
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channels = ["C2","C4";"C1","C3"]; %EEG channels to average. First row for left hand and second row for right hand
hands = ["L","R"];
colors = ["#a3a3a3","#3b3b3b"];

addpath(pdir)
sid = strcat('A',strtrim(string(num2str(sn(:))))); %This is hardcoded for our study. Change accordingly if necessary
filenames = strtrim(string(ls(strcat(pdir,'\',sid,'*',filesfx,'*.set'))));
ichans = zeros(length(channels),2);

figure('Units','inches','Position',[5,5,5,5],'Name',sid,'NumberTitle','off');
set(gcf,'Color','w')
mytile = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
for ihand = 1:length(hands)
    nexttile(mytile,ihand)
    EEG = pop_loadset('filename',char(filenames(ihand)),'filepath',pdir);
    istart = [1;EEG.ntrials1+1];
    iend = [EEG.ntrials1;EEG.ntrials2+EEG.ntrials1];
    %Find the channel indices
    for ichan = 1:size(channels,2)
        allchannels = struct2cell(EEG.chanlocs);
        chanlist = cellfun(@(x) any(strcmp(x,channels(ihand,ichan))),allchannels);
        ichans(ichan,ihand) = find(vertcat(chanlist(1,:,:)));
    end
    %Plot the data
    for isession = 1:2
        datamt = mean(EEG.data(:,:,istart(isession):iend(isession)),3);
        datamtc = mean(datamt(ichans(:,ihand),:),1);
        plot(EEG.times,-datamtc,'Color',colors(isession),'LineWidth',1.5)
        hold on
    end
    set(gca,'FontSize',16)
    title(strcat(hands(ihand),"H: ",strjoin(channels(ihand,:),', ')))
    xlim([10,300])
    xlabel('Time (ms)')
    ylabel('Voltage (\muV)')
    legend(["Pre","Post"],'location','southeast')
    legend boxoff
end