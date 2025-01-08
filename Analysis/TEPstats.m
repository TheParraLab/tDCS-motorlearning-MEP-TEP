function TEPstats(ihand,ptime,schans,type,S,chanlocs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Returns linear mixed effects stats for specified TEP data
%   
%   Input:
%   ihand -     hand index (L/R)
%   ptime -     peak time in ms post-TMS
%   schans -    EEG channels to analyze
%   type -      model type: 
%               'dNCSL' or 'dNCSR' for performance gain on left or right 
%               hand, respectively, or 
%               'dMEP' for post/pre MEP ratio
%   S -         structure array containing the data
%   chanlocs -  EEG channel info
%
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ng = 3;             %Number of groups
nsg = 40;           %Number of subjects per group
hands = ["L","R"];
intensity = [0;4;6]; %Stimulation intensity (mA)

nchan = length(schans);         %Number of EEG channels
prepostMEP = NaN(nsg*2,1,ng);   %Pre- and post-stimulation MEP amplitudes
dncs = NaN(nsg,ng);             %Performance gain
sn = NaN(nsg,ng);               %Subject IDs
dTEP = NaN(nsg,nchan,ng);       %TEP change
dMEP = NaN(nsg,ng);             %MEP change
Icol = reshape((intensity*ones(1,nsg))',nsg*ng,1);

for ig = 1:ng
    eval(strcat('prepostMEP(:,:,ig)=[S(ig).MEP',hands(ihand),'Pre;S(ig).MEP',hands(ihand),'Post];'))
    for j = 1:nchan
        allchannels = struct2cell(chanlocs);
        chanlist = cellfun(@(x) any(strcmp(x,schans(j))),allchannels);
        ichan = find(vertcat(chanlist(1,:,:)));
        eval(strcat('dTEP(:,j,ig)=diff(S(ig).P',num2str(ptime),'(:,ichan,ihand,:),1,4);'))
        eval(strcat('dMEP(:,ig)=S(ig).MPP',hands(ihand),';'))
    end
    sn(:,ig) = S(ig).sid;
    switch type
        case {'dNCSL','dNCSR'}
            switch type
                case 'dNCSL'
                    dncs(:,ig) = S(ig).ncs_diff(:,1); %Use task 1 for left hand
                case 'dNCSR'
                    dncs(:,ig) = S(ig).ncs_diff(:,3); %Use task 3 for right hand
            end            
    end
end

switch type
    case {'dNCSL','dNCSR'}
        dTEPcol = reshape(mean(dTEP,2,'omitnan'),[ng*nsg,1]);
        dncscol = reshape(dncs,[ng*nsg,1]);
        sncol = reshape(sn,[ng*nsg,1]);
        iremove = logical(sum([isnan(dTEPcol),isnan(dncscol)],2)); %Remove all rows with NaNs
        dTEPcol(iremove) = [];
        dncscol(iremove) = [];
        Icol(iremove) = [];
        sncol(iremove) = [];
        tbl = table(dTEPcol,dncscol,Icol,sncol,'VariableNames',{'dTEP','dNCS','Intensity','Subject'});
        tbl.Subject = categorical(tbl.Subject);
        lme = fitlme(tbl,'dTEP~Intensity*dNCS+(1|Subject)');
        disp(strcat("F(",num2str(lme.anova.DF1(4)),",",num2str(lme.anova.DF2(4)),") = ",num2str(lme.anova.FStat(4)),", p = ",num2str(lme.anova.pValue(4))))
    case 'dMEP'
        dTEPcol = reshape(mean(dTEP,2,'omitnan'),[ng*nsg,1]);
        dMEPcol = reshape(dMEP,[ng*nsg,1]);
        sncol = reshape(sn,[ng*nsg,1]);
        iremove = logical(sum([isnan(dTEPcol),isnan(dMEPcol)],2)); %Remove all rows with NaNs
        dTEPcol(iremove) = [];
        dMEPcol(iremove) = [];
        Icol(iremove) = [];
        sncol(iremove) = [];
        tbl = table(dTEPcol,dMEPcol,Icol,sncol,'VariableNames',{'dTEP','dMEP','Intensity','Subject'});
        tbl.Subject = categorical(tbl.Subject);
        lme = fitlme(tbl,'dTEP~Intensity*dMEP+(1|Subject)');
        disp(strcat("F(",num2str(lme.anova.DF1(4)),",",num2str(lme.anova.DF2(4)),") = ",num2str(lme.anova.FStat(4)),", p = ",num2str(lme.anova.pValue(4))))
end