%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script is used for quantitative analysis and visualization.
%   Figures are numbered according to the manuscript on PCI RR, available
%   here: https://osf.io/a42uy
% 
%   Copyright (c) 2024 Gavin Hsu and Zhenous Hadi Jafari
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear

maindir = input('Enter the path to your main directory: ',"s"); %Path to main directory
eeglabpath = input('Enter the path to your EEGLAB directory: ',"s"); %Path to EEGLAB

cd(maindir)
addpath bayesFactor/ %Place bayesFactor package here, under main directory
addpath Resources/
addeeglab(eeglabpath)

%Load data
load('combinedData.mat','S') %Preprocessed behavioral, MEP, and TEP data stored in a structure array
MEPdata = load('MEPdata.mat','data'); %MEP data, including traces
MEPdata = MEPdata.data;
MEPcell = struct2cell(MEPdata);
allfn = fieldnames(MEPdata);
TEPfile = load('TEPdata.mat','chanlocs','t','nchan'); %TEP parameters
nchan = TEPfile.nchan; %Number of channels in EEG
EEGt = TEPfile.t; %Time for EEG recording
chanlocs = TEPfile.chanlocs; %EEG channel locations
allchannels = struct2cell(chanlocs);

%Constants
ns = sum(cell2mat({S(:).n})); %Total number of subjects
ng = 3; %Number of groups
nsg = floor(ns/ng); %Number of subjects per group
pt = [33,45,55,100,175]'; %TEP peak times to evaluate
hands = ["L","R"]; %Left/Right hand
sessions = ["Pre","Post"]; %Pre/Post tDCS

%Colors and marker types for figures
color = [0.24,0,0.74;0.89,0.2,0.56;1,0.65,0]; %Foreground colors
color2 = ["#b999ff","#f3a5cd","#ffdb99"]; %Shade colors
markertypes = ["o","square","hexagram"]; %Markers for 0,4,6 mA, respectively

%Combine data across groups (for correlation stats)
allwpm = zeros(ng*nsg,1); %Typing speed
allsat = zeros(ng*nsg,1); %Saturated NCS in last nstrials trials
allear = zeros(ng*nsg,1); %Early NCS in first nstrials trials
allstart = zeros(ng*nsg,1); %First nonzero NCS
alldiff = zeros(ng*nsg,1); %Difference between early and saturated NCS
alldiffR = zeros(ng*nsg,1); %Difference between early and saturated NCS RH
allmncst = zeros(ng*nsg,1); %Mean NCS across trials
allmncstR = zeros(ng*nsg,1); %Mean NCS across trials RH
allmtst = zeros(ng*nsg,1); %Mean tapping speed across trials
alla_m = zeros(ng*nsg,1); %Mean accuracy across trials
allmep = NaN(ng*nsg,2); %Mean post/pre MEP ratio: subject, left/right
allmepPre = NaN(ng*nsg,1); 
allnkp = zeros(ng*nsg,1); %Total number of keypress throughout whole task
allr = zeros(ng*nsg,1); %Sensation ratings at beginning of stimulation
allTOD = zeros(ng*nsg,1); %Time of day
allage = zeros(ng*nsg,1);
allsex = cell(ng*nsg,1);
for ihand = 1:length(hands)
    for isession = 1:length(sessions)
        eval(strcat('allMEP',hands(ihand),sessions(isession),'=NaN(ng*nsg,1);')) %Mean MEP amplitude across trials
    end
end

for ig = 1:ng
    for ip = 1:length(pt)
        eval(strcat('allP',num2str(pt(ip)),'((ig-1)*nsg+1:ig*nsg,:,:,:)=S(ig).P',num2str(pt(ip)),';'))
    end
    allwpm((ig-1)*nsg+1:ig*nsg) = S(ig).wordspm;
    allsat((ig-1)*nsg+1:ig*nsg) = S(ig).ncs_sat_m(:,:,1);
    allear((ig-1)*nsg+1:ig*nsg) = S(ig).ncs_ear_m(:,:,1);
    allmncst((ig-1)*nsg+1:ig*nsg) = S(ig).mncst(:,1);
    allmncstR((ig-1)*nsg+1:ig*nsg) = S(ig).mncst(:,3);
    allstart((ig-1)*nsg+1:ig*nsg) = S(ig).ncs_start;
    alldiff((ig-1)*nsg+1:ig*nsg) = S(ig).ncs_diff(:,1);
    alldiffR((ig-1)*nsg+1:ig*nsg) = S(ig).ncs_diff(:,3);
    allmep((ig-1)*nsg+1:ig*nsg,1) = S(ig).MPPL;
    allmep((ig-1)*nsg+1:ig*nsg,2) = S(ig).MPPR;
    allmepPre((ig-1)*nsg+1:ig*nsg,1) = S(ig).MEPLPre;
    allmepPre((ig-1)*nsg+1:ig*nsg,2) = S(ig).MEPRPre;
    allr((ig-1)*nsg+1:ig*nsg) = S(ig).r(1,:)';
    allmtst((ig-1)*nsg+1:ig*nsg) = S(ig).mtst(:,1);
    alla_m((ig-1)*nsg+1:ig*nsg) = S(ig).a_m(:,:,1);
    allnkp((ig-1)*nsg+1:ig*nsg) = S(ig).nkp(1,:);
    allTOD((ig-1)*nsg+1:ig*nsg) = S(ig).TOD;
    allage((ig-1)*nsg+1:ig*nsg) = S(ig).age;
    allsex((ig-1)*nsg+1:ig*nsg) = S(ig).sex;
    for ihand = 1:length(hands)
        for isession = 1:length(sessions)
            eval(strcat('allMEP',hands(ihand),sessions(isession),'((ig-1)*nsg+1:ig*nsg)=S(ig).MEP',hands(ihand),sessions(isession),';'))
        end
    end
end

%% Statistics
%Copy data for analysis table
data = zeros(nsg,7,ng); %Array to be converted into a table for analysis
itask = 1; %Task number (1-4)
for ig = 1:ng
    data(:,1,ig) = S(ig).mncst(:,itask);
    data(:,2,ig) = S(ig).wordspm;
    data(:,3,ig) = S(ig).sid;
    data(:,4,ig) = S(ig).r(1,:)';
    data(:,5,ig) = S(ig).MPPL;
    data(:,6,ig) = S(ig).ncs_diff(:,itask);
    data(:,7,ig) = S(ig).ncs_start';
    data(:,8,ig) = squeeze(S(ig).ncs_sat_m(:,:,itask))';
end

alpha = 0.05; %Significance level
NCScol = reshape(data(:,1,:),[ng*nsg,1]); %NCS
TScol = reshape(data(:,2,:),[ng*nsg,1]); %Typing Speed
Rcol = reshape(data(:,4,:),[ng*nsg,1]); %Sensation Rating (beginning)
MPPLcol = reshape(data(:,5,:),[ng*nsg,1]); %MEP post/pre ratio left hand
Dcol = reshape(data(:,6,:),[ng*nsg,1]); %NCS gain
startcol = reshape(data(:,7,:),[ng*nsg,1]); %Trial 1 NCS
satcol = reshape(data(:,8,:),[ng*nsg,1]); %Saturated NCS
current = [zeros(nsg,1);4*ones(nsg,1);6*ones(nsg,1)]; %Current
subj = reshape(data(:,3,:),[ng*nsg,1]); %Subject ID
tbl_ncs = table(NCScol,TScol,Rcol,MPPLcol,Dcol,startcol,satcol,current,subj,'VariableNames',{'NCS','TS','Sensation','MEP','Delta','Start','Sat','Current','Subject'});
iremove = find(isnan(tbl_ncs.NCS));
tbl_ncs(iremove,:) = []; %Remove rows with missing NCS data

%Fit various linear models and calculate Bayes factors

%H2: MEP - Linear model with current as a continuous effect
lmMEP = fitlm(tbl_ncs,'MEP~Current');
bf10_MEP = bf.anova(tbl_ncs(~isnan(tbl_ncs.MEP),:),'MEP~Current');
bf01_MEP = 1/bf10_MEP; %Bayes factor for current intensity effect
tintercept = lmMEP.Coefficients.tStat(1);
pintercept = lmMEP.Coefficients.pValue(1);

%H3: NCS vs. MEP - Linear Mixed Effects model with MEP as a continuous effect
lmeNCSMEP = fitlme(tbl_ncs,'NCS~MEP+(1|Subject)');
bf10_H3 = bf.anova(tbl_ncs(~isnan(tbl_ncs.MEP),:),'NCS~MEP+(1|Subject)');
bf01_H3 = 1/bf10_H3;

%Exploratory: NCS vs. Sensation: Linear Mixed Effects model with sensation as a continuous effect
lme = fitlme(tbl_ncs,'NCS~Sensation*Current+(1|Subject)');
bf10_full = bf.anova(tbl_ncs,'NCS~Sensation*Current+(1|Subject)');
bf10_current = bf.anova(tbl_ncs,'NCS~Current+Sensation:Current+(1|Subject)');
bf10_sensation = bf10_full/bf10_current;
bf01_sensation = 1/bf10_sensation;

%Exploratory: NCS gain vs. MEP - Linear Mixed Effects model with MEP as a continuous effect
lmedNCSMEP = fitlme(tbl_ncs,'Delta~MEP+(1|Subject)');
bf10_dNCSMEP = bf.anova(tbl_ncs(~isnan(tbl_ncs.MEP),:),'Delta~MEP+(1|Subject)');
bf01_dNCSMEP = 1/bf10_dNCSMEP;

%Exploratory: Change in NCS
lmD = fitlm(tbl_ncs,'Delta~Current');

%Exploratory: Typing speed
lmTS = fitlm(tbl_ncs,'TS~Current');

%Exploratory: Baseline performance
lmStart = fitlm(tbl_ncs,'Start~Current');

%Exploratory: Saturated performance
lmSat = fitlm(tbl_ncs,'Sat~Current');

%H1: NCS - current as a continuous effect (Stage 1)
%Linear model with typing speed as covariate
lm = fitlm(tbl_ncs,'NCS~Current+TS'); %equivalent to: anovan(NCScol,{current,TScol},'display','off','continuous',[1,2])
p_cov = lm.anova.pValue(2);
bf10 = bf.anova(tbl_ncs,'NCS~Current'); %Bayes factor for current intensity effect

%Linear model without covariate
lm2 = fitlm(tbl_ncs,'NCS~Current'); %equivalent to: anovan(NCScol,current,'display','off','continuous',1);
p_nocov = lm2.anova.pValue(1);

%H1: NCS - current as a categorical effect (Stage 2)
if p_cov > alpha
    %Linear model with typing speed as covariate
    tbl_ncs.Current = categorical(tbl_ncs.Current); %Change current to a categorical effect
    lm3 = fitlm(tbl_ncs,'NCS~Current+TS'); %equivalent to: anovan(NCScol,{current,TScol},'display','off','continuous',2)
    p_cov2 = lm3.anova.pValue(2);
    
    %Linear model without covariate
    lm4 = fitlm(tbl_ncs,'NCS~Current'); %equivalent to: anovan(NCScol,current,'display','off','continuous',1);
    p_nocov2 = lm4.anova.pValue(1);

    %H1: NCS - pairwise comparison between 4 and 6mA (Stage 3)
    if p_cov2 < alpha/2
        [~,~,stats] = anovan(NCScol,{current,TScol},'display','off','continuous',2);
        mc = multcompare(stats,'Dimension',1,'Display','off'); %Post-hoc Tukey HSD
        pmc = mc(3,6); %Row 3 for groups 2 vs. 3 and Column 6 for p-value
        if pmc < alpha/3
            disp("H1: Reversing Effect")
        else
            bf10 = bf.ttest(data(:,1,2),data(:,1,3)); %Bayes factor for 4 vs. 6mA
            disp("Saturating Effect")
            disp(stracat("BF01 = ",num2str(1/bf10)))
        end
    else
        disp("H1: No Effect")
        disp(strcat("BF01 = ",num2str(1/bf10)))
    end
else
    disp("H1: Monotonic Effect")
end

%% Fit LME models for TEP Stats
% LME dTEP~Intensity*dNCS+(1|Subject)
TEPstats(1,45,["C2","C4"],'dNCSL',S,chanlocs)
TEPstats(1,55,["C2","C4"],'dNCSL',S,chanlocs)
TEPstats(2,45,["C1","C3"],'dNCSL',S,chanlocs)
TEPstats(2,55,["C1","C3"],'dNCSL',S,chanlocs)
% LME dTEP~Intensity*dMEP+(1|Subject)
TEPstats(2,45,["C1","C3"],'dMEP',S,chanlocs)
TEPstats(2,55,["C1","C3"],'dMEP',S,chanlocs)

%% Subject Demographics
SAS = readtable('demographics.csv');
meanage = mean(table2array(SAS(:,2)));
SDage = std(table2array(SAS(:,2)));
minage = min(table2array(SAS(:,2)));
maxage = max(table2array(SAS(:,2)));
nF = sum(strcmp(table2array(SAS(:,3)),'F'));
nM = sum(strcmp(table2array(SAS(:,3)),'M'));

%% Figure 3a: NCS vs. Trial
figure('Name','Figure 3a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
h = zeros(3,1);
hold on
xts = (1:36)'; %36 trials

for ig = 1:ng
    fill([xts;flipud(xts)],[S(ig).mncss(:,:,1)-S(ig).encss(:,:,1);flipud(S(ig).mncss(:,:,1)+S(ig).encss(:,:,1))],color(ig,:),'Linestyle','none','FaceAlpha',0.3);
    h(ig) = line(xts,S(ig).mncss(:,:,1),'Color',color(ig,:),'LineWidth',2);
end
xlabel('Trial')
ylabel('N Correct Sequences')
xlim([1,36])
legend(h,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','southeast')
legend box off
set(gcf,'color','w')
removeMargins

%% Figure 3b: Violin Plot NCS
figure('Name','Figure 3b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on

hm = zeros(3,1);
j = 1;
for ig = 1:ng
    xv = ig;
    yv = S(ig).mncst(:,j);
    cv = color(ig,:);
    mv = markertypes(ig);
    hm(ig) = vplot(xv,yv,cv,mv,1,50);
end

legend(hm,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','southoutside','FontSize',18,'NumColumns',3)
legend box off
xlim([0,4])
xticks(1:ng)
xticklabels([])
ylim([0,10])
ylabel('N Correct Sequences')
set(gcf,'color','w')
removeMargins

%% Figures 4 & S5: MEP Traces - Separated by Group + Violin Plot MEP
tMEP = cell(length(hands),length(sessions));
itMEP = cell(length(hands),length(sessions),size(MEPcell,3));
ifull = cell(length(hands),length(sessions),size(MEPcell,3));
linestyles = [":","-"];
hgroup = zeros(ng,length(sessions));
for ihand = 1:length(hands)
    for isession = 1:length(sessions)
        fnt = strcat('t',hands(ihand),'H',sessions(isession));
        ifnt = find(strcmp(allfn,fnt));
        alltMEP = MEPcell(ifnt,:);
        tMEP{ihand,isession} = unique(round(cat(1,alltMEP{:}),4));
        for isub = 1:size(MEPcell,3)-1
            [~,itMEP{ihand,isession,isub},ifull{ihand,isession,isub}] = intersect(round(alltMEP{isub},4),round(tMEP{ihand,isession},4));
        end
    end
end

tmin = 20;
tmax = 40;
pad = 7;
fignames = ["Figure 4","Figure S5"];
for ihand = 1:length(hands)
    figure('Name',fignames(ihand),'NumberTitle','off','Units','inches','Position',[5,5,9,4],'Color','white')
    tOut = tiledlayout(1,7,'Padding','tight','TileSpacing','compact');

    ticklocs = zeros(6,1);
    nvalid = zeros(ng,1);
    nexttile(tOut,[1,4])
    
    for ig = 1:ng
        for isession = 1:length(sessions)
            fntrace=strcat('traces',hands(ihand),'H',sessions(isession));
            ifntrace = find(strcmp(allfn,fntrace));
            alltraces = MEPcell(ifntrace,:);
    
            hold on
            MEPs = NaN(size(tMEP{ihand,isession},1),nsg);
            for isubject = 1:nsg
                if ~isnan(S(ig).iMEP(isubject)) && size(alltraces{S(ig).iMEP(isubject)},1)>1000 && size(alltraces{S(ig).iMEP(isubject)},2)>10
                    trace = median(alltraces{S(ig).iMEP(isubject)},2);
                    MEPs(ifull{ihand,isession,S(ig).iMEP(isubject)},isubject) = trace(itMEP{ihand,isession,S(ig).iMEP(isubject)});
                    tplot = 1000*tMEP{ihand,isession};
                    MEPplot = -MEPs(:,isubject);
                    itmin = find(tplot>=tmin,1);
                    itmax = find(tplot>tmax,1)-1;
                    hss = plot(tplot(itmin:itmax)+(pad+tmin)*(ig-1),MEPplot(itmin:itmax),linestyles(isession),'Color',color2(ig),'LineWidth',0.5);
                    nvalid(ig) = nvalid(ig) + 1;
                end
            end
            grouptrace = -mean(MEPs,2,'omitnan');
            hgroup(ig,isession) = plot(tplot(itmin:itmax)+(pad+tmin)*(ig-1),grouptrace(itmin:itmax),linestyles(isession),'Color',color(ig,:),'LineWidth',2);
            set(gca,'FontName','Arial','LineWidth',2,'FontSize',16)
        end
        for isession = 1:length(sessions)
            uistack(hgroup(ig,isession),'top')
        end
        ticklocs((ig-1)*2+1) = tplot(itmin)+(pad+tmin)*(ig-1);
        ticklocs((ig-1)*2+2) = tplot(itmax)+(pad+tmin)*(ig-1);
    end
    xlabel('Time (ms)','FontSize',20)
    ylabel('Voltage (\muV)','FontSize',20)
    ylim([-300,800])
    yticks([-300,0,800])
    xticks(ticklocs)
    xlim([ticklocs(1),ticklocs(end)])
    xticklabels([repmat([20,40],1,3)])
    
    nexttile(tOut,[1,3])
    hold on
    hm = zeros(3,1);
    nvalid = zeros(ng,1);
    for ig = 1:ng
        xv = ig;
        eval(strcat('yv = S(ig).MPP',hands(ihand),';'))
        cv = color(ig,:);
        mv = markertypes(ig);
        hm(ig) = vplot(xv,yv,cv,mv,1,50);
        eval(strcat('nvalid(ig)=sum(~isnan(S(ig).MPP',hands(ihand),'));'))
    end
    yline(1,'--')
    lgd = legend(hm,{['0mA n=',num2str(nvalid(1))],['4mA n=',num2str(nvalid(2))],['6mA n=',num2str(nvalid(3))]},'Location','southoutside','FontSize',16,'NumColumns',3);
    lgd.Box = 'off';

    ylim([0,4])
    xlim([0.5,3.5])
    xticks(1:ng)
    xticklabels([])
    ylabel('Post/Pre MEP Ratio')
    set(gca,'FontName','Arial','LineWidth',2,'FontSize',20)
end

%% Fig. 5a: NCS vs. MEP
figure('Name','Figure 5a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
nvalid = zeros(ng,1);
hm = zeros(3,1);
for ig = 1:ng
    hm(ig) = scatter(S(ig).MPPL,S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled');
    hold on
    nvalid(ig) = sum(~isnan(S(ig).MPPL));
end
lgd = legend(hm,{['0mA n=',num2str(nvalid(1))],['4mA n=',num2str(nvalid(2))],['6mA n=',num2str(nvalid(3))]},'Location','southeast','FontSize',18,'NumColumns',1,'EdgeColor','none','Color',[0.9,0.9,0.9]);
[r,p] = corr(allmep(:,1),allmncst,'Type','Pearson','rows','complete');
df = sum(~isnan(allmep(:,1))) - 2;
disp(strcat("r(",num2str(df),") = ",num2str(round(r,3,'significant')),", p = ",num2str(round(p,3,'significant'))))
ylabel('N Correct Sequences')
xlabel('Post/Pre MEP Ratio')
set(gcf,'color','w')
removeMargins

%% Fig. 5b: Delta NCS vs. MEP
figure('Name','Figure 5b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).MPPL,S(ig).ncs_diff(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allmep(:,1),alldiff,'Type','Pearson','rows','complete');
df = sum(~isnan(allmep(:,1)))-2;
disp(strcat("r(",num2str(df),") = ",num2str(round(r,3,'significant')),", p = ",num2str(round(p,3,'significant'))))
ylabel('\Delta NCS')
xlabel('Post/Pre MEP Ratio')
set(gcf,'color','w')
removeMargins

%% Figure 6: FTT Performance vs. Typing Speed
figure('Name','Figure 6a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).wordspm,S(ig).mncst(:,1),24,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(allmncst,[allwpm,ones(nsg*ng,1)]);
xfit = [min(allwpm),max(allwpm)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(allwpm,allmncst,'Type','Pearson');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('Typing Speed (WPM)')
set(gcf,'color','w')
removeMargins

figure('Name','Figure 6b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).wordspm,S(ig).ncs_start,24,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(allstart,[allwpm,ones(nsg*ng,1)]);
xfit = [min(allwpm),max(allwpm)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(allwpm,allstart,'Type','Pearson');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylim([0,5])
ylabel('Trial 1 NCS')
xlabel('Typing Speed (WPM)')
set(gcf,'color','w')
removeMargins

figure('Name','Figure 6c','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).wordspm,S(ig).ncs_diff(:,1),24,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(alldiff,[allwpm,ones(nsg*ng,1)]);
xfit = [min(allwpm),max(allwpm)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(allwpm,alldiff,'Type','Pearson');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('\Delta NCS')
xlabel('Typing Speed (WPM)')
set(gcf,'color','w')
removeMargins

%% Figure 7a: Sensation Ratings
figure('Name','Figure 7a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on
for ig = 1:ng
    errorbar(S(ig).mr,S(ig).er,'Color',color(ig,:),'LineWidth',2,'Marker',markertypes(ig),'MarkerSize',10,'MarkerFaceColor',color(ig,:))
end
legend(['0mA n=',num2str(size(S(1).r,2))],['4mA n=',num2str(size(S(2).r,2))],['6mA n=',num2str(size(S(3).r,2))],'Location','northeast','FontSize',20)

xlim([0.5,3.5])
ylim([0,10])
xticks(1:3)
xticklabels(["Beginning","Middle","After"])
xlabel('Stimulation Time Point')
ylabel('Sensation Rating')
legend box off

set(gcf,'color','w')
ax = gca;
set(ax,'FontName','Arial','LineWidth',2,'FontSize',20)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2) - 3*ti(4);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) + ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Figure 7b: Correlation between sensation and task performance
fig = figure('Name','Figure 7b','NumberTitle','off','Units','inches','Position',[5,5,5,4]);
gn = ["0mA";"4mA";"6mA"];

for ig = 1:ng
    scatter(S(ig).r(1,:),S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allmncst,allr,'Type','Pearson','rows','complete');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))

xlabel('Beginning Sensation Level','FontSize',20)
ylabel('Average NCS')

set(gcf,'color','w')
removeMargins

%% Figure 8: Scatter Plot NCS
figure('Name','Figure 8','NumberTitle','off','Units','inches','Position',[5,5,8,4])
hold on

d = 0.9; %displacement between box plots within group
xspread = 0.4;
h = zeros(3,1);
hm = zeros(3,1);
xpos = zeros(3,4);

for ig = 1:ng
    xpos(ig,:) = 1:4:16;
    xpos(ig,:) = xpos(ig,:)+(d*(ig-1));
    nvalid = S(ig).n-sum(isnan(S(ig).mncst),1);
    S(ig).encst = std(S(ig).mncst,0,1,'omitnan')./sqrt(nvalid); %SEM
    S(ig).cincs = S(ig).encst.*tinv(0.95,nvalid-1); %95 pct confidence interval
end

for j = 1:4 %task
    for ig = 1:ng %group
        h(ig) = plot(xpos(ig,j)+xspread/2,mean(S(ig).mncst(:,j),1,'omitnan'),'_','Color',color(ig,:),'MarkerSize',26,'LineWidth',2);
        xbox = [xpos(ig,j)-xspread/2;xpos(ig,j)+xspread*1.5];
        ybox = repmat(mean(S(ig).mncst(:,j),1,'omitnan'),2,1);
        fill([xbox;flipud(xbox)],[ybox-S(ig).encst(j);flipud(ybox+S(ig).encst(j))],color(ig,:),'Linestyle','none','FaceAlpha',0.3);
        xscat = rand(S(ig).n,1)*xspread;
        xscat = xscat + xpos(ig,j);
        hm(ig) = scatter(xscat,S(ig).mncst(:,j),24,color(ig,:),markertypes(ig),'filled');
    end
end

fill([4.5,4.5,17,17],[0,10,10,0],[0,0,0],'Linestyle','none','FaceAlpha',0.05);
legend(hm,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','northoutside','FontSize',18,'NumColumns',3)
legend box off
xlim([0,17])
xticks(2.5:4:16)
xticklabels(["L:S1","L:S1","R:S2","L:S3"])
ylabel('N Correct Sequences')
set(gcf,'color','w')
removeMargins

%% Figures 9 & S9: TEP Plots - averaged across sessions
st = 10; %starting time
et = 300; %ending time
it1 = find(EEGt==st); %index for starting time
it2 = find(EEGt==et); %index for ending time
hw = 5;
cvals = [-5,5];
%colormap for voltage
cmv = [linspace(0.039,1,128),linspace(1,0.76,128);linspace(0.039,1,128),linspace(1,0.039,128);linspace(0.76,1,128),linspace(1,0.039,128)]';
nvalid = zeros(ng,1);
channels = ["C2","C4";"C1","C3"];

fignames = ["Figure 9","Figure S9"];
for ihand = 1:length(hands)
    figure('Name',fignames(ihand),'NumberTitle','off','Units','inches','Position',[5,5,10,6],'NumberTitle','off');
    set(gcf,'Color','w')
    mytile = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
    nexttile(mytile,1)
    for ichan = 1:size(channels,2)
        chanlist = cellfun(@(x) any(strcmp(x,channels(ihand,ichan))),allchannels);
        ichans(ichan,ihand) = find(vertcat(chanlist(1,:,:)));
    end
    for ig = 1:ng   
        for ichan = 1:nchan
            eval(strcat('x=S(ig).TEP',hands(ihand),'(:,ichan,it1+1:it2,1);'))
            hchan = plot(EEGt(it1+1:it2),squeeze(mean(x,1,'omitnan')),'Color',color2(ig));
            uistack(hchan,'bottom')
            hold on
        end
        eval(strcat('x=S(ig).TEP',hands(ihand),'(:,ichans(:,ihand),it1+1:it2,1);'))
        nvalid(ig) = S(ig).n-sum(isnan(mean(x,[2,3])),1);
        x = mean(x,2,'omitnan');
        hac(ig) = plot(EEGt(it1+1:it2),squeeze(median(x,1,'omitnan')),'Color',color(ig,:),'LineWidth',2);
        set(gca,'FontSize',16,'LineWidth',2,'Box','off')
        xlim([st,et])
        ylim([-10,10])
    end
    for ipt = 1:length(pt)
        peakline = xline(pt(ipt),'--');
    end
    xlabel('Time (ms)','FontSize',16)
    ylabel('Voltage (\muV)','FontSize',16)
    legend(hac,{['0mA n=',num2str(nvalid(1))],['4mA n=',num2str(nvalid(2))],['6mA n=',num2str(nvalid(3))]},'Location','northeast','FontSize',16,'NumColumns',3,'Box','off')
    tile2 = tiledlayout(mytile,1,length(pt),'Padding','tight','TileSpacing','none');
    tile2.Layout.Tile = 2;
    for ipt = 1:length(pt)
        nexttile(tile2,ipt)
        iTEP = find(EEGt<pt(ipt)+hw&EEGt>pt(ipt)-hw);
        for ig = 1:ng
            eval(strcat('plotdata(ig,:) = squeeze(mean(S(ig).TEP',hands(ihand),'(:,:,iTEP,:),[1,3,4],''omitnan''));'))
        end
        topoplot(mean(plotdata,1),chanlocs,'whitebk','on');
        colormap(cmv)
        clim(cvals)
        title(strcat(num2str(pt(ipt))," ms"))
        set(gca,'FontSize',14)
        axis xy
    end
    cb = colorbar;
    cb.Layout.Tile = 'west';
    cb.Label.String = 'Voltage (\muV)';
    cb.Label.FontSize = 16;
    cb.Label.Position(1) = 3;
end

%% Figures 10 & S10: Change in TEP with Time Course
channels = ["C2","C4"];
ichans = zeros(size(channels,2),1);
for ichan = 1:size(channels,2)
    chanlist = cellfun(@(x) any(strcmp(x,channels(ichan))),allchannels);
    ichans(ichan) = find(vertcat(chanlist(1,:,:)));
end
for ip = 1:length(pt)
    eval(strcat('ptt2P',num2str(pt(ip)),'=zeros(nchan,2);')) %p-value for two-sample t-test between pre and post TEP
    eval(strcat('t2P',num2str(pt(ip)),'=zeros(nchan,2);')) %t-statistic for two-sample t-test
    for ichan = 1:nchan
        for ihand = 1:length(hands)
            eval(strcat('[psrP',num2str(pt(ip)),'(ichan,ihand),~,stats]=signrank(allP',num2str(pt(ip)),'(:,ichan,ihand,2),allP',num2str(pt(ip)),'(:,ichan,ihand,1));'))
            eval(strcat('zP',num2str(pt(ip)),'(ichan,ihand)=stats.zval;'))
        end
    end
end

st = 10; %starting time
et = 300; %ending time
it1 = find(EEGt==st); %index for starting time
it2 = find(EEGt==et); %index for ending time
hw = 5;
intensity = [0;4;6];
dTEP = NaN(nsg,size(channels,2),length(EEGt),ng);

fignames = ["Figure 10","Figure S10"];
for ihand = 1:length(hands)
    h = figure('Name',fignames(ihand),'NumberTitle','off','Units','inches','Position',[5,5,12,6],'NumberTitle','off');
    movegui('center')
    set(gcf,'Color','w')
    tOut = tiledlayout(h,2,1,'TileSpacing','tight','Padding','tight');
    tile1 = tiledlayout(tOut,1,1,'TileSpacing','none','Padding','tight');
    tile1.Layout.Tile = 1;
    ax1 = nexttile(tile1);
    hold on
    for ig = 1:ng
        eval(strcat('TEP=S(ig).TEP',hands(ihand),'(:,ichans,:,:);'))
        dTEP(:,:,:,ig) = diff(TEP,1,4);
        nvalid(ig) = S(ig).n-sum(isnan(mean(dTEP(:,:,:,ig),[2,3])),1);
        x = mean(dTEP(:,:,:,ig),2,'omitnan');
        edTEP = squeeze(iqr(x,1)/2);
        x = squeeze(median(x,1,'omitnan'));
        hac(ig) = plot(EEGt(it1+1:it2),x(it1+1:it2),'Color',color(ig,:),'LineWidth',2);
    end
    dTEPcol = reshape(mean(dTEP,[2,3]),nsg*ng,1);
    ymax = 5;
    ymin = -5;
    Icol = reshape((intensity*ones(1,nsg))',nsg*ng,1);
    iremove = isnan(dTEPcol);
    Icol(iremove) = [];
    tsig = [];

    yline(0,'--')
    for ipt = 1:length(pt)
        iTEP = find(EEGt<pt(ipt)+hw&EEGt>pt(ipt)-hw);
        dTEPcol = reshape(mean(dTEP(:,:,iTEP,:),[2,3]),nsg*ng,1);
        dTEPcol(iremove) = [];
        p = kruskalwallis(dTEPcol,Icol,'off');
        if p<0.05/length(pt)
            tsig = [tsig;pt(ipt)];
        end
    end
    if ~isempty(tsig)
        iborder = [0;find(diff(tsig)>1);length(tsig)];
        for iblock = 1:length(iborder)-1
            tshade = [tsig(iborder(iblock)+1)-hw;tsig(iborder(iblock+1))+hw;flipud([tsig(iborder(iblock)+1)-hw;tsig(iborder(iblock+1))+hw])];
            yshade = [ymin;ymin;ymax;ymax];
            hpatch = patch(tshade,yshade,'k','Linestyle','none','FaceAlpha',0.3);
            uistack(hpatch,'bottom')
        end
    end
    
    xlabel('Time (ms)','FontSize',16)
    ylabel('\Delta TEP (\muV)','FontSize',16)
    leg1 = legend(ax1,hac,{['0mA n=',num2str(nvalid(1))],['4mA n=',num2str(nvalid(2))],['6mA n=',num2str(nvalid(3))]},'Location','northeast','FontSize',16,'NumColumns',3,'Box','off');
    set(ax1,'FontName','Arial','FontSize',16,'LineWidth',2)
    xlim([st,et])
    ylim([ymin,ymax])
    %Add a second legend
    ax2 = axes('Parent',tile1,'Visible','off');
    leg2 = legend(ax2,hpatch,strcat('Kruskal-Wallis p<',num2str(0.05/length(pt))),'Location','southeast','FontSize',16,'Box','off');

    tile2 = tiledlayout(tOut,1,length(pt),'TileSpacing','none','Padding','tight');
    tile2.Layout.Tile = 2;
    for ip = 1:length(pt)
        nexttile(tile2,ip)
        eval(strcat('topoplot(zP',num2str(pt(ip)),'(:,ihand),chanlocs,''whitebk'',''on'',''emarker2'',{find(psrP',num2str(pt(ip)),'(:,ihand)<0.05),''*'',''m'',5});'))
        clim([-5,5])
        title(strcat(num2str(pt(ip))," ms"))
        set(gca,'FontSize',14)
        axis xy
    end
    cb1 = colorbar;
    cb1.Layout.Tile = 'west';
    cb1.Label.String = 'z-statistic';
    cb1.Label.FontSize = 16;
    cb1.Label.Position(1) = 3;
    title(tile2,'Wilcoxon Signed Rank Test Post-Pre TEP','FontSize',16,'FontWeight','bold')
end

%% Figure 11: TEP Topoplots for select peaks and measuress
prepostMEP = NaN(ng*nsg*2,2);
for ihand = 1:length(hands)
    eval(strcat('prepostMEP(:,ihand)=[allMEP',hands(ihand),'Pre;allMEP',hands(ihand),'Post];')) %Pre AND Post MEPs
end
for ip = 1:length(pt)
    eval(strcat('prepostTEP',num2str(pt(ip)),'=NaN(ng*nsg*2,2);'))
    eval(strcat('rP',num2str(pt(ip)),'=zeros(nchan,2);')) %Pearson correlation coefficient for post-pre difference vs. post/pre MEP ratio
    eval(strcat('pP',num2str(pt(ip)),'=zeros(nchan,2);')) %p-value for Pearson correlation
    eval(strcat('ptt2P',num2str(pt(ip)),'=zeros(nchan,2);')) %p-value for two-sample t-test between pre and post TEP
    eval(strcat('t2P',num2str(pt(ip)),'=zeros(nchan,2);')) %t-statistic for two-sample t-test
    for ichan = 1:nchan
        for ihand = 1:length(hands)
            eval(strcat('prepostTEP',num2str(pt(ip)),'(:,ihand)=[allP',num2str(pt(ip)),'(:,ichan,ihand,1);allP',num2str(pt(ip)),'(:,ichan,ihand,2)];')) %Pre AND Post TEPs
            eval(strcat('dTEP=diff(allP',num2str(pt(ip)),'(:,ichan,ihand,:),1,4);')) %Post-Pre difference
            eval(strcat('[rP',num2str(pt(ip)),'(ichan,ihand),pP',num2str(pt(ip)),'(ichan,ihand)]=corr(dTEP,allmep(:,ihand),''Type'',''Pearson'',''rows'',''complete'');'))
            eval(strcat('[~,ptt2P',num2str(pt(ip)),'(ichan,ihand),~,stats]=ttest(allP',num2str(pt(ip)),'(:,ichan,ihand,1),allP',num2str(pt(ip)),'(:,ichan,ihand,2));'))
            eval(strcat('t2P',num2str(pt(ip)),'(ichan,ihand)=stats.tstat;'))
            eval(strcat('[rPrePostP',num2str(pt(ip)),'(ichan,ihand),pPrePostP',num2str(pt(ip)),'(ichan,ihand)]=corr(prepostMEP(:,ihand),prepostTEP',num2str(pt(ip)),'(:,ihand),''Type'',''Pearson'',''rows'',''complete'');'))
            eval(strcat('[rNCSP',num2str(pt(ip)),'(ichan,ihand),pNCSP',num2str(pt(ip)),'(ichan,ihand)]=corr(alldiff,dTEP,''Type'',''Pearson'',''rows'',''complete'');'))
            eval(strcat('[rNCSRP',num2str(pt(ip)),'(ichan,ihand),pNCSRP',num2str(pt(ip)),'(ichan,ihand)]=corr(alldiffR,dTEP,''Type'',''Pearson'',''rows'',''complete'');'))
        end
    end
end

figure('Name','Figure 11','NumberTitle','off','Units','inches','Position',[5,5,12,6],'NumberTitle','off');
movegui('center')
set(gcf,'Color','w')
colormap('jet')
tOut = tiledlayout(1,5,'TileSpacing','none','Padding','tight');
tIn1 = tiledlayout(tOut,2,2,'TileSpacing','none','Padding','tight');
tIn1.Layout.Tile = 1;
tIn1.Layout.TileSpan = [1,2];
ihand = 1;
for ip = 2:3
    nexttile(tIn1,(ihand-1)*2+ip-1)
    eval(strcat('topoplot(rNCSP',num2str(pt(ip)),'(:,ihand),chanlocs,''whitebk'',''on'',''emarker2'',{find(pNCSP',num2str(pt(ip)),'(:,ihand)<0.05),''*'',''m'',5});'))
    clim([-0.5,0.5])
    title(strcat(hands(ihand),"H t = ",num2str(pt(ip))," ms"))
    set(gca,'FontSize',14)
    axis xy
end
ihand = 2;
for ip = 2:3
        nexttile(tIn1,(ihand-1)*2+ip-1)
        eval(strcat('topoplot(rNCSRP',num2str(pt(ip)),'(:,ihand),chanlocs,''whitebk'',''on'',''emarker2'',{find(pNCSRP',num2str(pt(ip)),'(:,ihand)<0.05),''*'',''m'',5});'))
        clim([-0.5,0.5])
        title(strcat(hands(ihand),"H t = ",num2str(pt(ip))," ms"))
        set(gca,'FontSize',14)
        axis xy
end
title(tIn1,'\DeltaTEP vs. \DeltaNCS','FontSize',16,'FontWeight','bold')

tIn2 = tiledlayout(tOut,2,2,'TileSpacing','none','Padding','tight');
tIn2.Layout.Tile = 3;
tIn2.Layout.TileSpan = [1,2];
for ihand = 1:length(hands)
    for ip = 2:3
        nexttile(tIn2,(ihand-1)*2+ip-1)
        eval(strcat('topoplot(rP',num2str(pt(ip)),'(:,ihand),chanlocs,''whitebk'',''on'',''emarker2'',{find(pP',num2str(pt(ip)),'(:,ihand)<0.05),''*'',''m'',5});'))
        clim([-0.5,0.5])
        title(strcat(hands(ihand),"H t = ",num2str(pt(ip))," ms"))
        set(gca,'FontSize',14)
        axis xy
    end
end
title(tIn2,'\DeltaTEP vs. MEP Ratio','FontSize',16,'FontWeight','bold')

tIn3 = tiledlayout(tOut,2,1,'TileSpacing','none','Padding','tight');
tIn3.Layout.Tile = 5;
tIn3.Layout.TileSpan = [1,1];
for ihand = 1:length(hands)
    for ip = 4
        nexttile(tIn3,(ihand-1)+ip-3)
        eval(strcat('topoplot(rPrePostP',num2str(pt(ip)),'(:,ihand),chanlocs,''whitebk'',''on'',''emarker2'',{find(pPrePostP',num2str(pt(ip)),'(:,ihand)<0.05),''*'',''m'',5});'))
        clim([-0.5,0.5])
        title(strcat(hands(ihand),"H t = ",num2str(pt(ip))," ms"))
        set(gca,'FontSize',14)
        axis xy
    end
end
title(tIn3,'TEP vs. MEP','FontSize',16,'FontWeight','bold')

cb1 = colorbar;
cb1.Layout.Tile = 'east';
cb1.Label.String = 'Pearson Correlation Coefficient r';
cb1.Label.FontSize = 16;
cb1.Label.FontWeight = 'bold';

%% Figure S1: NCS vs. Speed, NCS vs. Accuracy, MEP vs. NKP, MEP vs. Typing Speed
figure('Name','Figure S1a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).mtst(:,1),S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(allmncst,[allmtst,ones(nsg*ng,1)]);
xfit = [min(allmtst),max(allmtst)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(allmtst,allmncst,'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('Average Speed (kp/s)')
set(gcf,'color','w')
removeMargins

figure('Name','Figure S1b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).a_m(:,:,1),S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(allmncst,[alla_m,ones(nsg*ng,1)]);
xfit = [min(alla_m),max(alla_m)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(alla_m,allmncst,'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('Average Accuracy')
set(gcf,'color','w')
removeMargins

figure('Name','Figure S1c','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).a_m(:,:,1),S(ig).mtst(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(alla_m,allmtst,'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Average Speed (kp/s)')
xlabel('Average Accuracy')
set(gcf,'color','w')
removeMargins

figure('Name','Figure S1d','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).ncs_diff(:,1),S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled')
    hold on
end
b = regress(allmncst,[alldiff,ones(nsg*ng,1)]);
xfit = [min(alldiff),max(alldiff)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(alldiff,allmncst,'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('\Delta NCS')
set(gcf,'color','w')
removeMargins

%% Figure S2a-b: Time of Day
figure('Name','Figure S2a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on

wpmAM = allwpm(~allTOD);
wpmPM = allwpm(allTOD==1);
yy = {wpmAM,wpmPM};
for iTOD = 1:2
    y = yy{iTOD};
    y = y(~isnan(y));
    [f1,xi1] = ksdensity(y);
    f1 = f1/max(f1)*0.4;
    patch([-f1,fliplr(f1)]+iTOD-1,[xi1,fliplr(xi1)],'k','FaceAlpha',0.2,'LineStyle','none')
    plot(iTOD-1,mean(y),'_','Color','k','MarkerSize',100,'LineWidth',2);
    for ig = 1:ng
        iiTOD = find(S(ig).TOD==iTOD-1);
        C = arrayfun(@(arr) arr-xi1,S(ig).wordspm(iiTOD),'UniformOutput',false);
        C = cell2mat(C);
        [~,iw] = min(abs(C),[],2);
        w = fliplr(f1)+f1;
        xscat = (rand(length(iiTOD),1)-0.5).*w(iw)'*0.5;
        scatter(S(ig).TOD(iiTOD)+xscat,S(ig).wordspm(iiTOD),30,color(ig,:),markertypes(ig),'filled')
    end
end
[~,p,~,stats] = ttest2(wpmAM,wpmPM);
disp(strcat("t(",num2str(stats.df),") = ",num2str(stats.tstat),", p = ",num2str(p)))
ylabel('Typing Speed (WPM)')
xlabel('Time')
set(gcf,'color','w')
xlim([-0.5,1.5])
xticks([0,1])
xticklabels(["Morning","Afternoon"])
removeMargins

figure('Name','Figure S2b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on
mncstAM = allmncst(~allTOD);
mncstPM = allmncst(allTOD==1);
yy = {mncstAM,mncstPM};
for iTOD = 1:2
    y = yy{iTOD};
    y = y(~isnan(y));
    [f1,xi1] = ksdensity(y);
    f1 = f1/max(f1)*0.4;
    patch([-f1,fliplr(f1)]+iTOD-1,[xi1,fliplr(xi1)],'k','FaceAlpha',0.2,'LineStyle','none')
    plot(iTOD-1,mean(y),'_','Color','k','MarkerSize',100,'LineWidth',2);
    for ig = 1:ng
        iiTOD = find(S(ig).TOD==iTOD-1);
        C = arrayfun(@(arr) arr-xi1,S(ig).mncst(iiTOD),'UniformOutput',false);
        C = cell2mat(C);
        [~,iw] = min(abs(C),[],2);
        w = fliplr(f1)+f1;
        xscat = (rand(length(iiTOD),1)-0.5).*w(iw)'*0.5;
        scatter(S(ig).TOD(iiTOD)+xscat,S(ig).mncst(iiTOD),30,color(ig,:),markertypes(ig),'filled')
    end
end
[~,p,~,stats] = ttest2(mncstAM,mncstPM);
disp(strcat("t(",num2str(stats.df),") = ",num2str(stats.tstat),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('Time')
set(gcf,'color','w')
xlim([-0.5,1.5])
xticks([0,1])
xticklabels(["Morning","Afternoon"])
removeMargins

%% Figure S2c-d: Subject Sex
figure('Name','Figure S2c','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on

wpmM = allwpm(strcmp(allsex,'M'));
wpmF = allwpm(strcmp(allsex,'F'));
yy = {wpmM,wpmF};
sexes = ['M','F'];
for i = 1:2
    y = yy{i};
    y = y(~isnan(y));
    [f1,xi1] = ksdensity(y);
    f1 = f1/max(f1)*0.4;
    patch([-f1,fliplr(f1)]+i-1,[xi1,fliplr(xi1)],'k','FaceAlpha',0.2,'LineStyle','none')
    plot(i-1,mean(y),'_','Color','k','MarkerSize',100,'LineWidth',2);
    for ig = 1:ng
        ii = find(strcmp(S(ig).sex,sexes(i)));
        C = arrayfun(@(arr) arr-xi1,S(ig).wordspm(ii),'UniformOutput',false);
        C = cell2mat(C);
        [~,iw] = min(abs(C),[],2);
        w = fliplr(f1)+f1;
        xscat = (rand(length(ii),1)-0.5).*w(iw)'*0.5;
        scatter(i-1+xscat,S(ig).wordspm(ii),30,color(ig,:),markertypes(ig),'filled')
    end
end
[~,p,~,stats] = ttest2(wpmM,wpmF);
disp(strcat("t(",num2str(stats.df),") = ",num2str(stats.tstat),", p = ",num2str(p)))
ylabel('Typing Speed (WPM)')
xlabel('Sex')
set(gcf,'color','w')
xlim([-0.5,1.5])
xticks([0,1])
xticklabels(["Male","Female"])
removeMargins

figure('Name','Figure S2d','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on
mncstM = allmncst(strcmp(allsex,'M'));
mncstF = allmncst(strcmp(allsex,'F'));
yy = {mncstM,mncstF};
sexes = ['M','F'];
for i = 1:2
    y = yy{i};
    y = y(~isnan(y));
    [f1,xi1] = ksdensity(y);
    f1 = f1/max(f1)*0.4;
    patch([-f1,fliplr(f1)]+i-1,[xi1,fliplr(xi1)],'k','FaceAlpha',0.2,'LineStyle','none')
    plot(i-1,mean(y),'_','Color','k','MarkerSize',100,'LineWidth',2);
    for ig = 1:ng
        ii = find(strcmp(S(ig).sex,sexes(i)));
        C = arrayfun(@(arr) arr-xi1,S(ig).mncst(ii)','UniformOutput',false);
        C = cell2mat(C);
        [~,iw] = min(abs(C),[],2);
        w = fliplr(f1)+f1;
        xscat = (rand(length(ii),1)-0.5).*w(iw)'*0.5;
        scatter(i-1+xscat,S(ig).mncst(ii),30,color(ig,:),markertypes(ig),'filled')
    end
end
[~,p,~,stats] = ttest2(mncstM,mncstF);
disp(strcat("t(",num2str(stats.df),") = ",num2str(stats.tstat),", p = ",num2str(p)))
ylabel('Average NCS')
xlabel('Sex')
set(gcf,'color','w')
xlim([-0.5,1.5])
xticks([0,1])
xticklabels(["Male","Female"])
removeMargins

%% Figure S2e-f: Subject Age

figure('Name','Figure S2e','NumberTitle','off','Units','inches','Position',[5,5,5,4])
nvalid = zeros(ng,1);
hm = zeros(3,1);
for ig = 1:ng
    hm(ig) = scatter(S(ig).age,S(ig).wordspm(:,1),30,color(ig,:),markertypes(ig),'filled');
    hold on
    nvalid(ig) = sum(~isnan(S(ig).age));
end
b = regress(allwpm,[allage,ones(nsg*ng,1)]);
xfit = [min(allage),max(allage)];
yfit = [b(2)+b(1)*xfit(1),b(2)+b(1)*xfit(2)];
plot(xfit,yfit,'--','Color','k','LineWidth',1.5)
[r,p] = corr(allage(:,1),allwpm,'Type','Pearson','rows','complete');
df = sum(~isnan(allage(:,1))) - 2;
disp(strcat("r(",num2str(df),") = ",num2str(round(r,3,'significant')),", p = ",num2str(round(p,3,'significant'))))
ylabel('Typing Speed (WPM)')
xlabel('Age (years)')
set(gcf,'color','w')
removeMargins

figure('Name','Figure S2f','NumberTitle','off','Units','inches','Position',[5,5,5,4])
nvalid = zeros(ng,1);
hm = zeros(3,1);
for ig = 1:ng
    hm(ig) = scatter(S(ig).age,S(ig).mncst(:,1),30,color(ig,:),markertypes(ig),'filled');
    hold on
    nvalid(ig) = sum(~isnan(S(ig).age));
end
[r,p] = corr(allage(:,1),allmncst,'Type','Pearson','rows','complete');
df = sum(~isnan(allage(:,1))) - 2;
disp(strcat("r(",num2str(df),") = ",num2str(round(r,3,'significant')),", p = ",num2str(round(p,3,'significant'))))
ylabel('Average NCS')
xlabel('Age (years)')
set(gcf,'color','w')
removeMargins

%% Figure S3: Violin Plot Delta NCS
figure('Name','Figure S3','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on

hm = zeros(3,1);
j = 1;
for ig = 1:ng %group
    xv = ig;
    yv = S(ig).ncs_diff(:,j);
    cv = color(ig,:);
    mv = markertypes(ig);
    hm(ig) = vplot(xv,yv,cv,mv,1,50);
end

legend(hm,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','southoutside','FontSize',18,'NumColumns',3)
legend box off
xlim([0,4])
xticks(1:ng)
xticklabels([])
ylim([0,7])
ylabel('\Delta NCS')
set(gcf,'color','w')
removeMargins

%% Figure S4: MEP vs. Typing
figure('Name','Figure S4a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).nkp(1,:),S(ig).MPPL,30,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allnkp,allmep(:,1),'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Post/Pre MEP Ratio')
xlabel('Number of Key Presses')
set(gcf,'color','w')
removeMargins

figure('Name','Figure S4b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).wordspm,S(ig).MEPLPre,30,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allwpm,allmepPre(:,1),'Type','Pearson','rows','complete');
title(strcat("r = ",num2str(r),", p = ",num2str(p)))
ylabel('Pre-tDCS MEP (\muV)')
xlabel('Typing Speed (WPM)')
set(gcf,'color','w')
removeMargins

%% Figure S6a: Violin Plot Saturated NCS
figure('Name','Figure S6a','NumberTitle','off','Units','inches','Position',[5,5,5,4])
hold on

hm = zeros(3,1);
j = 1;
for ig = 1:ng %group
    xv = ig;
    yv = squeeze(S(ig).ncs_sat_m(:,:,j))';
    cv = color(ig,:);
    mv = markertypes(ig);
    hm(ig) = vplot(xv,yv,cv,mv,1,50);
end

legend(hm,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','southoutside','FontSize',18,'NumColumns',3)
legend box off
xlim([0,4])
xticks(1:ng)
xticklabels([])
ylim([0,10])
ylabel('Saturated NCS')
set(gcf,'color','w')
removeMargins

%% Figure S6b: MEP vs. Trial 1 NCS
figure('Name','Figure S6b','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).MPPL,S(ig).ncs_start',24,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allmep(:,1),allstart,'Type','Pearson','rows','complete');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))
% title('Saturated NCS vs. Typing Speed')
ylabel('Trial 1 NCS')
xlabel('Post/Pre MEP Ratio')
set(gcf,'color','w')
removeMargins

%% Figure S6c: MEP vs. Saturated NCS
figure('Name','Figure S6c','NumberTitle','off','Units','inches','Position',[5,5,5,4])
for ig = 1:ng
    scatter(S(ig).MPPL,S(ig).ncs_sat_m(:,:,1),24,color(ig,:),markertypes(ig),'filled')
    hold on
end
[r,p] = corr(allmep(:,1),allsat,'Type','Pearson','rows','complete');
disp(strcat("r = ",num2str(r),", p = ",num2str(p)))
% title('Saturated NCS vs. Typing Speed')
ylabel('Saturated NCS')
xlabel('Post/Pre MEP Ratio')
set(gcf,'color','w')
removeMargins

%% Figure S7: Sensation questionnaire
severityData = readtable('Resources/severityrows.xlsx');

% tDCS symptoms and sensation
symptoms = {'Headache', 'Neck Pain', 'Scalp Pain', 'Tingling', 'Burning Sensation', 'Skin Redness', 'Sleepiness', 'Trouble Concentrating', 'Acute Mood Changes'};

figure('Name','Figure S7','NumberTitle','off','Units', 'inches', 'Position', [4,0,8,7], 'Color', 'w'); 
mytile = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for i = 1:9
    nexttile(mytile,i) 
    hold on
    hm = zeros(ng, 1);
    
    for ig = 1:ng
        xv = ig;
        rowIndices = S(ig).sid; 
        yv = severityData{rowIndices, i+2}; 
        cv = color(ig, :);
        mv = markertypes(ig);
        hm(ig) = vplot(xv, yv, cv, mv,2,40);
    end

    % Aesthetics 
    xlim([0.5, 3.5]);
    ylim([1, 5]);
    yticks(1:5);
    ylabel(symptoms{i});
    xticks([])
    set(gca, 'FontName', 'Arial', 'LineWidth', 2, 'FontSize', 14); 

    if i == 8
        legend(hm, {'0mA n=40', '4mA n=40', '6mA n=40'}, 'Location', 'southoutside', 'FontSize', 16, 'NumColumns', 3);
        legend box off
    end
end

%% Figure S8: Scatter Plot Tapping Speed
figure('Name','Figure S8','NumberTitle','off','Units','inches','Position',[5,5,10,5])
hold on

d = 0.9; %displacement between box plots within group
xspread = 0.4;
h = zeros(3,1);
hm = zeros(3,1);
xpos = zeros(3,4);

for ig = 1:ng
    xpos(ig,:) = 1:4:16;
    xpos(ig,:) = xpos(ig,:)+(d*(ig-1));
    nvalid = S(ig).n-sum(isnan(S(ig).mtst),1);
    S(ig).etst = std(S(ig).mtst,0,1,'omitnan')./sqrt(nvalid); %SEM
end

for j = 1:4 %task
    for ig = 1:ng %group
        h(ig) = plot(xpos(ig,j)+xspread/2,mean(S(ig).mtst(:,j),1,'omitnan'),'_','Color',color(ig,:),'MarkerSize',26,'LineWidth',2);
        xbox = [xpos(ig,j)-xspread/2;xpos(ig,j)+xspread*1.5];
        ybox = repmat(mean(S(ig).mtst(:,j),1,'omitnan'),2,1);
        fill([xbox;flipud(xbox)],[ybox-S(ig).etst(j);flipud(ybox+S(ig).etst(j))],color(ig,:),'Linestyle','none','FaceAlpha',0.3);
        xscat = rand(S(ig).n,1)*xspread;
        xscat = xscat + xpos(ig,j);
        hm(ig) = scatter(xscat,S(ig).mtst(:,j),24,color(ig,:),markertypes(ig),'filled');
    end
end

fill([4.5,4.5,17,17],[0,8,8,0],[0,0,0],'Linestyle','none','FaceAlpha',0.05);
legend(hm,{['0mA n=',num2str(S(1).n)],['4mA n=',num2str(S(2).n)],['6mA n=',num2str(S(3).n)]},'Location','northoutside','FontSize',18,'NumColumns',3)
legend box off
xlim([0,17])
xticks(2.5:4:16)
xticklabels(["Concurrent Stim","Lasting Effect","Other Hand","Other Sequence"])
ylabel('Speed (kp/s)')
set(gcf,'color','w')
removeMargins