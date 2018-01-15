%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters for calcium channels and IP3 transient
% Global variables are accessed by Paper_Ca_ODE_TH file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
clear;

global gamma delta v_ip3r v_leak v_in k_out v_pmca k_pmca
global d1 d2 d3 d5 a2
global v_serca k_serca
global k_soc v_soc
global stim_time d_rise d_decay
global r_rise Amp

%% Pump/Receptor/Leak Parameters
gamma=5.4054;
 
% Leak for ER
v_leak=0.002; 

% Leak for Extracellular Space
v_in=0.05; 
k_out=1.2;

% IP3R Parameters
v_ip3r=0.222;
% Li-Rinzel Parameters
d1=0.13; d2=1.049; d3=943.4e-3; d5=0.08234; 
a2=0.04; %adjusted Li-Rinzel Parameter

% PMCA Terms
v_pmca=10; k_pmca=2.5;

% SOCC Terms
v_soc=1.57;k_soc=90;

% SERCA Terms
v_serca=0.9;k_serca=0.1;
 
% Sneyd Parameter
delta=0.2;

%% IP3 Parameters (reproduces a SP response) 
stim_time=20;
Amp=0.375; 
d_rise=17;
d_decay=97; 
r_rise=0.002; 

 