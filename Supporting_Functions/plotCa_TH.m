function [TrueTroughVal, TrueTroughLoc, PeakVal,...
    PeakLoc]=plotCa_TH(x, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns everything from the "PeaksAndTroughs" function as well as plots
%the x trace (varargin 1 is the time vector and is not required; also 
%varargin 2 is the TitleName and is not required)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>=2 && strcmp(varargin(2),'time')==1
    t=varargin{1};
else
    t=1:length(x);
end

plot(t, x,'b', 'linewidth', 2); xlim([0 t(end)]);
if nargin>=2 && strcmp(varargin(2),'title')==1
    TitleName=varargin(1);
    title(TitleName)
end
hold on

[TrueTroughVal, TrueTroughLoc, PeakVal,...
    PeakLoc] = PeaksAndTroughs_TH(x, 0.26, 0.03); 

%If we want to plot the peaks and troughs detected in the
%"PeaksAndTroughs_TH" function:
%NewPeakLoc=t(PeakLoc); NewTroughLoc=t(TrueTroughLoc);
%plot(NewPeakLoc,PeakVal,'xm','markersize',10)
%plot(NewTroughLoc,TrueTroughVal,'ok', 'markersize', 7)

