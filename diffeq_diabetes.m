function dxdt = diffeq_diabetes(t,x,p_internal,c,ti,xi,input)
%diffeq_diabetes - Differential equations of the E-DES model.
%
%   dxdt = diffeq_diabetes(t,x,p_internal,c,input)
%
%Inputs
% t:    current time point
% x:    vector with state variables (at current time t)
% p_internal: vector of model parameters
% c:    struct of model constants
%input: struct with model inputs (time of the meal, carbohydrate content of
%       the meal, timing and dose of exogenous insulin 
%
%Output
% dxdt: vector with derivatives of state variables

%History
% 5-Mar-2016 Check on integrator time window (Natal van Riel, TU/e)
% 5-Mar-2016 Help info added (Natal van Riel, TU/e)

% -- state variables
Mg    = x(1); %Mg: glucose mass in gut compartment 1[mg]
Gpl   = x(2); %Gpl: plasma glucose concentration [mmol/L]
Ipl   = x(3); %Ipl: plasma insulin concentration [mU/L]
Gint  = x(4); %Gint: integrated plasma glucose increase (int (Gpl)) [mmol/L]
Id1   = x(5);

% -- evaluating splined AAs % -- might be better alternatives %
Apl = interp1(ti, xi(:,1), t,'pchip');
Aplb = xi(1,1);
dApl = interp1(ti, xi(:,2),t,'pchip');

% -- model constants
f      = c.f;       %f [mmol/mg], conversion factor from mg/dl to mmol/l
vg     = c.vg;      %vg [L/kg], glucose distribution volume
gbliv  = c.gbliv;   %gbliv [mmol/L], basal liver glucose output
beta   = c.beta;    %beta [(mmol/L)/(microU/mL)], conversion from microU/ml to mmol/l
taui   = c.taui;    %taui [min], integration time constant
taud   = c.taud;    %taud [min], differentation time constant
Gthpl  = c.Gthpl;   %Gthpl [mmol/L], threshold for renal glucose removal
c1     = c.c1;      %c1 [1/min](previously k8); constant term in gren
t_integralwindow = c.t_integralwindow; %Lower bound of moving time window of Gint

% -- model input
D       = input.D; %D (total amount of carbohydrates ingested) [mg] 
t_meal_start = input.t_meal_start;% starting time of the meal [min, counted from 0 = 0:00]
Mb      = input.Mb; %Mb (body weight) [kg]

% -- model parameters
k1    = p_internal(1); %[1/min]
k2    = p_internal(2); %[1/min]
k3    = p_internal(3); %[1/min]
k4    = p_internal(4); %[1/min]
k5    = p_internal(5); %[1/min]
k6    = p_internal(6); %[1/min]
k7    = p_internal(7); %[1/min]
k8    = p_internal(8); %[1/min]
k9    = p_internal(9);  %[1/min]
k10   = p_internal(10); %[1/min]
k11   = p_internal(11); %[1/min]
k12   = p_internal(12); %[1/min]
k13   = p_internal(13); %[1/min]
sigma = p_internal(14); %[-]
KM    = p_internal(15); %[mmol/l]
Gb    = p_internal(16); %[mmol/l]
Iplb  = p_internal(17); %[microU/ml]


%% ------------ ODEs --------------------------------------------------------
%% dMg/dt -- glucose in gut 
if t>t_meal_start
    mgmeal  = sigma .* k1.^sigma .* (t-t_meal_start).^(sigma-1) .* exp(-(k1.*(t-t_meal_start)).^sigma) .* D ; 
else
    mgmeal = 0;
end

mgpl1    = k2.*Mg;
dMg_dt  = mgmeal - mgpl1;

%% dGpl/dt -- glucose in plasma
ggut    = k2 .*(f/(vg*Mb)) .* Mg; 
gliv    = gbliv - k3.*(Gpl-Gb) - k4.*beta.*(Id1) + k11.*(Apl-Aplb);

gnonit  = gbliv.*((KM+Gb)/Gb) .* (Gpl./(KM+Gpl));
git     = k5.*beta.*Id1.*(Gpl./(KM+Gpl)); 

if Gpl > Gthpl
    gren  = c1./(vg*Mb) .* (Gpl - Gthpl); 
else
    gren = 0;
end

dGpl_dt = gliv + ggut - gnonit - git - gren;

%% dIpl/dt -- insulin in plasma
% ipnc
global t_saved Gpl_saved

t_lowerbound = t - t_integralwindow;
if (t > t_integralwindow) && (length(t_saved)>1) && (length(t_saved) == length(Gpl_saved))
    Gpl_lowerbound = interp1(t_saved,Gpl_saved,t_lowerbound, 'spline');
else
    Gpl_lowerbound = Gpl_saved(1);  % is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
end
%Gpl_lowerbound = interp1(t_saved,Gpl_saved,t_lowerbound, 'spline');

dGint_dt = (Gpl-Gb) - (Gpl_lowerbound-Gb);
ipnc    = (beta.^-1).*(k6.*(Gpl-Gb) + (k7/taui).*Gint + (k7/taui).*Gb + (k8.*taud).*dGpl_dt) + k12.*(dApl) + k13.*(Apl-Aplb);

%iliv & iif
iliv    = k7.*(Gb./(beta*taui*Iplb)).*Ipl;
iif     = k9.*(Ipl-Iplb); 

dIpl_dt = ipnc  - iliv - iif;

%% Insulin interstitial fluid 
dId1_dt = iif - k10.*Id1;

%% -- catch an error where the timestep of the integration becomes too small
% MINSTEP = 1e-10; %Minimum step
% 
% persistent tprev
% 
% if isempty(tprev)
%     tprev = -inf;
% end
% timestep = t - tprev;
% tprev = t;
% 
% if (timestep > 0) && (timestep < MINSTEP)
%     error(['Stopped. Time step is too small: ' num2str(timestep)])
% end


%% -- output differential equations
dxdt = [dMg_dt
        dGpl_dt
        dIpl_dt
        dGint_dt
        dId1_dt];