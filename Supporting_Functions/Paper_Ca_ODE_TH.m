%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full ODE model with IP3R, SERCA, leak, ECS leak, PMCA out, and SOC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = Paper_Ca_ODE_TH(t, x0)

CaCyt = x0(1); % calcium concentration in the cytosol
CaT = x0(2);   % total free calcium concentration in the cytosol
h = x0(3);     % deactivating variable for the IP3R

%% Global Parameters
global gamma delta v_ip3r v_leak v_in k_out v_pmca k_pmca
global d1 d2 d3 d5 a2
global v_serca k_serca
global k_soc v_soc
global stim_time d_rise d_decay
global r_rise Amp

% finds the IP3 value at the current time
[ ip ] = ip_function_TH(d_rise, d_decay, r_rise, Amp, stim_time, t);

% Uses conservation to find the calcium in the ER
CaER = (CaT-CaCyt)*gamma;

%% Terms for Calcium Dynamics on ER:

% terms for IP3R
minf = ip/(ip + d1);
ninf = CaCyt/(CaCyt + d5);
j_ip3r = v_ip3r*minf^3*ninf^3*h^3*(CaER - CaCyt);

% terms for the h ODE
q2 = d2*(ip+d1)/(ip+d3);
tauh = 1/(a2*(q2+CaCyt));
hinf = q2/(q2+CaCyt);

% Leak Term
j_leak = v_leak*(CaER-CaCyt);

%SERCA Pump
%note Hill coefficient of 1.75 as opppose to 2
j_serca = v_serca*CaCyt^1.75/(CaCyt^1.75 + k_serca^1.75);


%% Terms for Calcium Dynamics on Plasma Membrane:

%PMCA pump
j_pmca = v_pmca*CaCyt^2/(k_pmca^2 + CaCyt^2);

%SOCC 
j_soc=v_soc*k_soc^4./(k_soc^4+CaER.^4);

%Leak Terms
j_out = k_out*CaCyt;
j_in = v_in;


%% ODEs
xdot(1) = j_ip3r-j_serca+j_leak + (j_in-j_out+j_soc-j_pmca)*delta; %ODE for [Ca]cyt
xdot(2) = (j_in-j_out+j_soc-j_pmca)*delta; %ODE for [Ca]ER
xdot(3) = (hinf - h)/tauh; %ODE for h

xdot=xdot';
end

