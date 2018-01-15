%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example file that simulates a calcium transients with a 
% specified IP3 input and collects transient characteristics/type
% Uses default channel parameters, which are found in Official_Params file
% Included IP3 parameters reproduce the simulations found in Fig. 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
addpath('./Supporting_Functions');

Official_Params_TH; % loads the default parameters

% uncomment the line of code to produce wanted response
% reproduces calcium transients in Fig. 2
% JCNS, DOI: 10.1007/s10827-017-0640-1
% Order of parameters: amp d_rise r_rise d_decay 179
ip3_params=[0.2, 10, 0.2, 97 ]; % single-peak
%ip3_params=[0.26 41 0.12 200]; % multi-peak 
%ip3_params=[0.375, 34, 0.002, 138]; % plateau
%ip3_params=[0.55, 39, 0.002, 179]; % long-lasting

% overwrites the IP3 global variables to the chosen values
Amp = ip3_params(1);
d_rise = ip3_params(2);
r_rise = ip3_params(3);
d_decay = ip3_params(4);

% sets the simulation time 
dt = 0.01;
sim_t = [0:dt:120];

% finds the IP3 transient
IP3 = ip_function_TH(d_rise, d_decay, r_rise, Amp, stim_time, sim_t);

% the simulation begins at steady state (based on default parameter values)
init_cond=[0.0865415 36.49084 0.6255124];

% runs the simulation using ode45
options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.1);
[sim_t,ca_transient] = ode15s(@Paper_Ca_ODE_TH, sim_t, init_cond, options);
% saves the calcium transient

CaCyt = ca_transient(:,1);
CaT =  ca_transient(:,2);
h =  ca_transient(:,3);

%% Plots the calcium transient and outputs response type 
[TrueTroughVal, TrueTroughLoc, PeakVal,PeakLoc]=plotCa_TH(CaCyt,...
    sim_t, 'time');
set(gca,'fontsize',16)
axis([0 max(sim_t) 0 max(CaCyt)])
xlabel('Time (sec)','fontsize',16)
ylabel('[Ca] (\muM)','fontsize',16)

%Prints the calcium response type
[Result, CaDur, CaAmount, CaLatency, StartOfResp,EndOfResp]=...
    FourCaResponseTypes_TH(TrueTroughVal,TrueTroughLoc, PeakVal,...
    PeakLoc, CaCyt, sim_t, 20);

%% Data collection and plot

% Stores characteristics about the IP3 transient
[IP3Amount, IP3PkLoc, IP3Amp, IP3TrueDur] = IP3Dyn_model_TH(IP3,sim_t);

% Stores characteristics about the Calcium transient
CaIP3PeakLatency=PeakLatency_TH(CaCyt, IP3);

CaAmp = max(PeakVal);

if strcmp(Result, 'NR')==0
    [CaRiseTime, CaDecTime] = CaRiseAndDecay_TH(CaCyt, dt);
else
    CaRiseTime=0; CaDecTime=0; CaAmp=0;
end

%Collects and saves information regarding the calcium and IP3 transients 
AllResults1=[d_decay, d_rise, r_rise, IP3Amp, IP3Amount,IP3TrueDur,...
    CaAmp, CaAmount, CaDur, CaLatency, CaRiseTime, CaDecTime];
AllResults2=[Result, CaIP3PeakLatency];

clearvars -except AllResults1 AllResults2

