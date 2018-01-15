function [IP3Amount, IP3PkLoc, IP3Amp, IP3TrueDur]=IP3Dyn_model_TH(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given the IP3 trace and corresponding time, it finds several IP3
%characteristics: total amount (area under curve), peak time, peak
%amplitude, and total IP3 duration (from start to where IP3 reached <0.005).
%Multiple things will need to change in this function if the baseline [IP3] 
%is not set at 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[IP3Amp, IP3PkLocTemp] = findpeaks(x);
IP3PkLoc = t(IP3PkLocTemp);

IP3AmountTemp = cumsum(x);
IP3Amount = IP3AmountTemp(end)/length(x)*t(end); %This is fine for IP3 amount
%because the baseline is 0. If that changes, then may need to change this (like Ca).
%%Or:
%IP3Amount = trapz(x)/length(x)*t(end);

if isempty(IP3PkLocTemp)==0
    t_end = t(IP3PkLocTemp + find(x(IP3PkLocTemp:end)<0.005, 1,'first') - 1);
else
    t_end =[];
end


t_start = t(find(x>0, 1,'first')); %1st point where IP3>0
if isempty(t_end)==0 %if there is an end to the IP3 dynamics
    IP3TrueDur = t_end - t_start;
else
    IP3TrueDur = t(end); %assume the IP3 is for the duration of the entire
    %simulation (including the time before any IP3 was applied!)
end

end
