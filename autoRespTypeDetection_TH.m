%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Marsa Taheri and Gregory Handy, 2016
% This code was used to simulate the mathematical model of Astrocyte 
% IP3-dependent Ca responses in 2 papers submitted in Nov 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code systematically steps through various IP3 traces and finds the 
% resulting Ca response dynamics. By entering a set of IP3 parameters, it 
% generates details about the resulting IP3 and Ca trace. As an example, 6
% total sets of IP3 parameters are chosen here (2 different Amp values, 1
% d_rise and d_decay value, and 3 r_rise values) for illustration. To
% generate the result for all 600 IP3 parameters (as discussed in Taheri et
% al.), use the values in the comments for these 4 IP3 parameters (Amp,
% d_rise, d_decay, r_rise).
%
% NOTE: The "pause" command on line 132 allows the user to view each Ca and
% IP3 trace before the code generates the next one (the user must press a
% key to unpause the code). This line may be commented out to quickly 
% generate the final results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
clear, clc, close all;

addpath('./Supporting_Functions');

Official_Params_TH; %Official_Params_TH;

dt = 0.001;
exp_t = [0:dt:290]'; %270 seconds after the input at t=20 s (defined by ip_input)

%% I.C.s:
%I.C. for the official params starting from steady state:
x0(1)=0.086541496665789; x0(2)=36.490839775010841;
x0(3)=0.625512446053023;

%% Run simulations:

ip_input=20; %IP3 input time

options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.1);

countAll = 0;
count_Amp=0;
for Amp = [0.55]; %or use linspace(0.2,0.9,5) for all Amp values
    count_drise = 0;
    count_Amp=count_Amp+1;
    
    LengthDrise = 5;%keep at 5 for all drise value, related to next line (i.e. length of d_rise values)
    for d_rise = 11; %or use linspace(1,41,LengthDrise) for all d_rise values
        count_drise=count_drise+1;

        if d_rise<8; r_rise_Vals = [0.002, 12];
        elseif d_rise<15; r_rise_Vals = [0.002, 0.44, 1.6];
        elseif d_rise<30; r_rise_Vals = [0.002, 0.12, 0.3, 1];
        elseif d_rise<40; r_rise_Vals = [0.002, 0.07, 0.15, 0.3, 0.8];
        else r_rise_Vals=[0.002, 0.04, 0.09, 0.15, 0.3, 0.8];
        end

        N_r_rise=length(r_rise_Vals);
        
        count_rrise = 0;
        for r_rise = r_rise_Vals; %or use "r_rise_Vals" to go through all possible r_rise values from above
            clear t x_sim CaCyt
            count_rrise = count_rrise + 1;
            
            s_0 = Amp/(1-exp(-r_rise*(d_rise)));
            
            d_decay_values = 179; %or use linspace(15,220,6) for all d_decay values
            for j=1:length(d_decay_values)
                d_decay=d_decay_values(j);
                countAll = countAll+1;

                
                [t,x_sim] = ode45(@Paper_Ca_ODE_TH,exp_t,x0,options);                
                
                ip_value = ip_function_TH(d_rise, d_decay, r_rise, Amp, ip_input, t);
                
%                 CaCyt_all(countAll, :) = x_sim(:,1); %used if need to save the traces
%                 IP3Cyt_all(countAll, :) = ip_value; %used if need to save the traces

                clf(figure(1)) %In order to show the last Ca response only; may comment out
                figure(1)
                plot(t, ip_value,'g', 'linewidth', 2)
                hold on
                
                %%%Finds IP3 characteristics:
                
                [IP3Amount, IP3PkLoc, IP3Amp, IP3TrueDur] = IP3Dyn_model_TH(ip_value,t);
                
                
                
                %%%The next several lines find Calcium characteristics:
                
                [TrueTroughVal, TrueTroughLoc, PeakVal,...
                    PeakLoc]=plotCa_TH(x_sim(:,1), t, 'time'); %returns the location 
                %indices; The troughs and peaks it finds are used below in
                %"FourCaResponseTypes_TH" to detect the Ca response type.
                
                CaIP3PeakLatency=PeakLatency_TH(x_sim(:,1), ip_value); %Latency of Ca peak
                %WRT IP3 peak, in terms of index number. To convert to time: t(CaIP3PeakLatency{1})
                
                CaAmp = max(PeakVal); %Calcium peak
                
                [Result, CaDur, CaAmount, CaInpLatency, StartOfResp,...
                    EndOfResp]=FourCaResponseTypes_TH(TrueTroughVal,...
                    TrueTroughLoc, PeakVal, PeakLoc, x_sim(:,1), t, ip_input); % Returns 
                %Calcium response type as well as its total duration, total amount 
                %(area under the curve), it's start and end points, and the latency of 
                %when the Ca response began WRT when the IP3 input began.
                
                if strcmp(Result, 'NR')==0
                    [CaRiseTime, CaDecTime] = CaRiseAndDecay_TH(x_sim(:,1), dt); 
                    %Finds the Ca rise/decay times
                else
                    CaRiseTime=0; CaDecTime=0; CaAmp=0; 
                    %so that they're not empty for the AllResults1 variable
                end                
                
                
                title({['Resp. Type= ',Result]; ['A= ', num2str(Amp),...
                    ', d_r= ', num2str(d_rise) ', r_r= ', num2str(r_rise),...
                    ', d_d= ', num2str(d_decay)]},...
                    'fontsize', 10)
                
                legend('IP3','Ca2+')
                xlabel('Time (sec)', 'fontsize', 10)
                ylabel('Conc (\muM)', 'fontsize', 10)
                set(gca,'fontsize', 10)
                
                
                %Use when saving data:
                AllResults1(countAll,:)=[d_decay, d_rise, r_rise, IP3Amp, IP3Amount,...
                    IP3TrueDur, CaAmp, CaAmount, CaDur, CaInpLatency, CaRiseTime, CaDecTime];
                AllResults2(countAll,:)=[Result, CaIP3PeakLatency];

                display('Press any key to unpause'); pause
                figure(2); plotCaAmountVsIP3Amount_TH %Example figure to generate

            end     
        end
    end    
end

%To save the data:
%save('DefaultParams_AllIP3', 'AllResults1','AllResults2')

figure(2);
set(gca,'fontsize', 11)
xlabel('Total IP3 Amount (\muM)', 'fontsize', 12)
ylabel('Total Ca2+ Amount (\muM)', 'fontsize', 12)

GenerateLegend_TH %Generate the Legend
