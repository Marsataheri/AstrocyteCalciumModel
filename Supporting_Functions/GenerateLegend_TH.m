%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script generates the figure legend for the 2nd figure in code 
%"autoRespTypeDetection_TH.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Order of Colors: SP, MP, Plat, LL
AllColors = [0.05,0.6,0.9; 0.7,0.1,0.1; 0,0,0; 0.8,0.8,0.9];
AllBorders = [0.15 0.4 1; 0.4 0 0; 0,0,0; 0.65,0.6,0.7];
AllFontSizes= 24; AllMarkerSizes = 22;

figure; hold on
plot(1,2.5,'s', 'markersize', AllMarkerSizes,'MarkerEdgeColor', AllBorders(1,:),...
        'MarkerFaceColor',AllColors(1,:))
text(1.1,2.5,'Single-Peak', 'fontsize', AllFontSizes)
plot(1,1.5,'s', 'markersize', AllMarkerSizes,'MarkerEdgeColor', AllBorders(2,:),...
        'MarkerFaceColor',AllColors(2,:))
text(1.1,1.5,'Multi-Peak', 'fontsize', AllFontSizes)
plot(1,2,'sk', 'markersize', AllMarkerSizes,...
        'MarkerFaceColor',AllColors(3,:))
text(1.1,2,'Plateau', 'fontsize', AllFontSizes)
plot(1,1,'s', 'MarkerEdgeColor', AllBorders(4,:),...
        'markersize', AllMarkerSizes,'MarkerFaceColor', AllColors(4,:))
text(1.1,1,'Long-Lasting', 'fontsize', AllFontSizes)


xlim([0.8, 1.8]); ylim([0.2, 3.3]);
set(gca, 'YTickLabel', [], 'YTick', [], 'XTickLabel', [],...
    'XTick', [], 'fontsize', 16)
set(gca, 'Visible','off')
