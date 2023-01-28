function [x0,input] = load_x0_input(sim,data,k,timespan)
%load_x0_input - Information to initialize a simulation.
%
%   [x0,input] = load_x0_input(sim,data,k,timespan)
%
%inputs:
% sim:      struct with information of a previous simulation
% fields    sim.Mg       glucose mass in gut [mg]
%           sim.Gpl      plasma glucose concentration [mmol/L]
%           sim.Ipl      plasma insulin concentration [mU/L]
%           sim.Uisc1    insulin concentration in intracellular compartment 1 [mU/L]
%           sim.Uisc2    insulin concentration in intracellular compartment 2 [mU/L]
%           sim.Gint     integrated plasma glucose increase (int (Gpl)) [mmol/L]
%           sim.Cpep     cpeptide mass in the body [nmol]
%           sim.Y        cpeptide mass in the interstitial fluid [nmol]
%           sim.Ipl_measured  measured insulin concentration (differs from Ipl because not all injected insulins are picked up by the assay [mU/L]
% data:     struct with data of one subject
% k:        indicates which meal {1,2,3}
% timespan: time span
%outputs:
% x0:       initial conditions (vector)
% input:    struct with information on the meal and insulin used

%History
% 5-Mar-2016 Help info added (Natal van Riel, TU/e)

%% -- set initial values for state variables equal to ending value of previous peak   
x0 = [sim.Mg(end)       %Mg: glucose mass in gut [mg]
      sim.Gpl(end)      %Gpl: plasma glucose concentration [mmol/L]
      sim.Ipl(end)      %Ipl: plasma insulin concentration [mU/L]
      sim.Gint(end)
      sim.Id1(end)
      ];     %Gint: integrated plasma glucose increase (int (Gpl)) [mmol/L]];   % measured insulin concentration (differs from Ipl because not all injected insulins are picked up by the assay [mU/L] 

%% -- set input: meal variables
input = struct();
input.D    = data.Dmeal{1,k};   %D (total amount of carbohydrates ingested) [mg]
input.t_meal_start = data.t_meal_start(1,k); % starting time of the meal [min, counted from 0 = 0:00]
input.Mb   = data.Mb;           %Mb (body weight) [kg]

end

  
