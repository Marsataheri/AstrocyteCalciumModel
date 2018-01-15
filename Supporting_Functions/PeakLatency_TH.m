function CaIP3PeakLatency=PeakLatency_TH(x,input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the index number of the difference between when the IP3 and Ca
%peaks occurred. Also, it describes which peak occurred first.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MaxCaPeak = find(x==max(x));
MaxInput = max(input);
AfterOrBefore = find(input==MaxInput) - MaxCaPeak;
if AfterOrBefore>0
    CaIP3PeakLatency = {AfterOrBefore, 'IP3 peak is LATER'};
elseif  AfterOrBefore<0
    CaIP3PeakLatency = {AfterOrBefore, 'IP3 peak is EARLIER'};
else
    CaIP3PeakLatency = {AfterOrBefore, 'Same time as IP3 peak'};
end
