%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Detects and returns step components with a slope with and order of
%   magnitude greater than -3
%
%   Inputs:
%   EEG -       EEGLAB EEG structure output from the fastICA step
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oom = findstepIC(EEG)
icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
icaact = double(icaact);
icaact = reshape(icaact, size(icaact,1), EEG.pnts, EEG.trials);
micaact = mean(icaact,3);
slope = zeros(EEG.nbchan,1);
x = 1:EEG.pnts;
x = [x',ones(EEG.pnts,1)];
for i = 1:EEG.nbchan
    b = regress(micaact(i,:)',x);
    slope(i) = b(1);
end
oom = find(log10(abs(slope))>-3); %find components with slope at an order of magnitude greater than -3
disp(strcat("Step components found: ",strjoin(strcat('IC',string(oom)),', ')))