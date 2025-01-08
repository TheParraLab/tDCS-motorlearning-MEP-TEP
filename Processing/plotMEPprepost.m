function plotMEPprepost(sn,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Visualizes MEPs for a given subject, comparing pre/post MEPs on 
%   different hands
%
%   Inputs:
%   sn -        Subject ID number, e.g. 33
%   data -      Structure array from MEPdata.mat
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check for subject data
sns = vertcat(data.sn);
isub = find(sns==sn,1);
if isempty(isub)
    disp("Subject data not found!")
else
    sessions = ["Pre","Post"];
    hands = ["LH","RH"];
    colors = ["#0072BD","#D95319"];
    figure('Units','inches','Position',[5,5,10,5]);
    set(gcf,'Color','w')
    mytile = tiledlayout(2,3,'TileSpacing','tight','Padding','compact');

    for ihand = 1:2
        nexttile(mytile,(ihand-1)*3+1)
        set(gca,'FontSize',16)
        ylabel('Voltage (\muV)')
        xlim manual
        xlim([0,0.05])
        title(strcat("Subject ",num2str(sns(isub))," ",hands(ihand)))
        hold on
        for isession = 1:2
            eval(strcat('t',num2str(isession),'=data(',num2str(isub),').t',hands(ihand),sessions(isession),';'))
            eval(strcat('avgtrace=mean(data(',num2str(isub),').traces',hands(ihand),sessions(isession),',2);'))
            eval(strcat('plot(t',num2str(isession),',-avgtrace,''Color'',colors(isession),''LineWidth'',1.5)'))
        end
        legend(sessions,'Location','northeast')
        legend box off
        nexttile(mytile,(ihand-1)*3+2)
        eval(strcat('imagesc(t1,1:size(data(',num2str(isub),').traces',hands(ihand),sessions(1),',2),-data(',num2str(isub),').traces',hands(ihand),sessions(1),''')'))
        xlim([0,0.05])

        set(gca,'FontSize',16)
        ylabel('Trial')
        title('Pre')
        nexttile(mytile,(ihand-1)*3+3)
        eval(strcat('imagesc(t2,1:size(data(',num2str(isub),').traces',hands(ihand),sessions(1),',2),-data(',num2str(isub),').traces',hands(ihand),sessions(2),''')'))
        xlim([0,0.05])
        set(gca,'FontSize',16)
        title('Post')
        xlabel(mytile,'Time (s)','FontSize',16)
    end
end