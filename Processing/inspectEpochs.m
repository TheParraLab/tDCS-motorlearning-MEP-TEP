function inspectEpochs(session,channel,edir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inspects epoched data for a specific channel, showing one epoch at a
%   time
%
%   Inputs:
%   session -   epoched recording session identifier, e.g. 'A33LPre'
%   channel -   the channel to view
%   edir -      path to directory with epoched EEG data
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = pop_loadset(strcat(session,'.set'),edir);
allchannels = struct2cell(EEG.chanlocs);
ichan = cellfun(@(x) any(strcmp(x,channel)),allchannels);
ichan = find(vertcat(ichan(1,:,:)));
figure('Units','inches','Position',[5,5,7,5]);
set(gcf,'Color','w')

for itrial = 1:size(EEG.data,3)
    plot(EEG.times,EEG.data(ichan,:,itrial),'LineWidth',1.5)
    title(strcat(channel," Trial ",num2str(itrial)))
    set(gca,'FontSize',16)
    xlabel('Time (ms)','FontSize',16)
    ylabel('Voltage (\muV)','FontSize',16)
    % xlim([-10,200])
    ylim([-100,100])
    pause
end