function [p,c]=load_pHealthy_c(data,glucose,insulin,i)

%load_pHealthy_c - Reference parameters for healthy subject model
%
%[p,c]=load_pHealthy_c()
%
%Inputs:
% data: struct with information on insulin used. Use [] for a normal simulation
% of a healthy subject
%Outputs:
% p: parameters (vector)
% c: constants (struct)

%History
% 5-Mar-2016 Help info added (Natal van Riel, TU/e)

% -- healthy parameter values from the best fitting model
p = [0.01    %k1
    0.28     %k2
    6.07e-3  %k3   
    2.35e-4  %k4 
    9.49e-2  %k5 
    1.93e-1  %k6 
    1.15     %k7  
    7.27     %k8 
    3.83e-2  %k9 
    2.84e-1  %k10
    5.0e-6   %k11 % AA parameter
    0.01	 %k12
	1.0e-4   %kl3 % AA parameter
    1.34     %sigma /index 14
    13.2     %Km /index 15
    glucose(i)%Gb /index 16
    insulin(i)]; %Ib/ index 17

%%
% -- set constants
c        = struct();
c.f      = 0.005551;    %f [mmol/mg], must be equal to mgdL_to_mmolL /10 in load_data.m
c.vg     = 17/70;       %vg [L/kg]
c.gbliv  = 0.043;       %gbliv [mmol/L]
c.beta   = 1;           %beta [(mmol/L)/(microU/mL)]
c.taui   = 31;          %taui [min]
c.taud   = 3;           %taud [min]
c.vi     = 13/70;       %distribution volume of insulin per kg bodymass, vi [L/kg]
c.Gthpl  = 9;           %Gthpl [mmol/L]
c.t_integralwindow = 30;%Lower bound of moving time window of Gint
c.c1     = 0.1;         %c1 [1/min](previously k8)

