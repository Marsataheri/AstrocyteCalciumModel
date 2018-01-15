function [CaRiseTime, CaDecTime, slope_up, slope_down] = CaRiseAndDecay_TH(CaTrace, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouputs the 10% to 90% rise and decay times and slopes of the calcium trace;
%dt is the time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline=CaTrace(2);

AboveBaseline=baseline*1.4;

[amp_max,time_of_max] = max(CaTrace);

first_trough = find(CaTrace>AboveBaseline,1,'first');

height = amp_max-baseline;
starting_amp = CaTrace(first_trough);
new_amp=height-starting_amp;

%Finds the 10% and 90% mark going up:
ninety_up = find(CaTrace > 0.9*new_amp + starting_amp, 1);
ten_up = find(CaTrace > 0.10*new_amp + starting_amp, 1);

%Finds the 10% and 90% mark going down:
ninety_down = find(CaTrace(time_of_max:end) > 0.9*new_amp + starting_amp, 1, 'last') + time_of_max;
ten_down = find(CaTrace(time_of_max:end) > 0.10*new_amp + starting_amp, 1, 'last') + time_of_max;

%Calculates the slope up and down:
slope_up = (CaTrace(ninety_up)-CaTrace(ten_up))...
    /(ninety_up-ten_up);
slope_down = (CaTrace(ninety_down)-CaTrace(ten_down))...
    /(ninety_down-ten_down);

%Calculates the rise time and decay time:
CaRiseTime = (ninety_up-ten_up)*dt;
CaDecTime = (ten_down-ninety_down)*dt;

end
