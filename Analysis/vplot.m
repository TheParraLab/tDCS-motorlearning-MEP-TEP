function h = vplot(x,y,color,marker,xspread,markersize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Creates a single violin plot and returns 'h' as a handle for the 
%   scatter plot.
%
%   Inputs:
%   x -         position along the horizontal axis to place the plot at
%   y -         data values along the vertical axis
%   color -     RGB triplet for the color of the violin
%   marker -    marker type for the scatter points
%   xspread -   double for controlling point spread width along the 
%               horizontal axis
%   markersize -size of the markers
% 
%   Copyright (c) 2024 Gavin Hsu
%   Parra Lab, CCNY, New York, NY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vw = 0.4;                   %Width of violin shape
y = y(~isnan(y));
[f1,xi1] = ksdensity(y);    %Kernel density estimate for y
f1 = f1/max(f1)*vw;         %Normalize violin width
C = arrayfun(@(arr) arr-xi1,y,'UniformOutput',false);
C = cell2mat(C);
[~,iw] = min(abs(C),[],2);  %Find points in xi1 closest to each value in y 
xscat = (rand(1,length(y))-0.5).*f1(iw)*xspread;    %Randomly scatter points within violin shape
patch([-f1,fliplr(f1)]+x,[xi1,fliplr(xi1)],color,'FaceAlpha',0.3,'LineStyle','none')
plot(x,mean(y),'_','Color',color,'MarkerSize',markersize,'LineWidth',2);
h = scatter(xscat+x,y,30,color,marker,'filled');