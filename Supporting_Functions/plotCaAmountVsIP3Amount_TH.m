%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the total Ca amount (CaAmount) as a function of total IP3
% amount (IP3Amount) with each symbol/square color showing the Ca response 
% type (SP, MP, Plat, or LL).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors3=[0.05 0.6 0.9; 0, 0, 0; 0.7, 0.1, 0.1; 0.8,0.8,0.9]; 
SizeOfMarker = 9;

hold on
if strcmp(Result,'SP')==1
    plot(IP3Amount,CaAmount,'s', 'MarkerEdgeColor', [0.15 0.4 1],...
        'markersize', SizeOfMarker, 'MarkerFaceColor',colors3(1,:)) 

elseif strcmp(Result,'MP')==1
    plot(IP3Amount,CaAmount,'s', 'MarkerEdgeColor', [0.4 0 0],...
    'markersize', SizeOfMarker, 'MarkerFaceColor',colors3(3,:))

elseif strcmp(Result,'Plat')==1
    plot(IP3Amount,CaAmount,'sk', 'markersize', SizeOfMarker,...
        'MarkerFaceColor',colors3(2,:))

elseif strcmp(Result,'LL')==1
plot(IP3Amount,CaAmount,'s', 'MarkerEdgeColor', [0.65,0.6,0.7],...
    'markersize', SizeOfMarker,'MarkerFaceColor', colors3(4,:))
end
