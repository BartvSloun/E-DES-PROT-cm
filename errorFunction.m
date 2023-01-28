function E = errorFunction(p0, pRest, namP,timespan,datGlu,datIns,timespan_i,datA_val,x0,c,input,sim) % removed pDat

for ind= 1:length(p0)
pRest(namP{ind})= p0(ind);
end

ptemp=pRest;
global t_saved Gpl_saved
t_saved = 0;
Gpl_saved = sim.Gpl;

ODE_model    = @diffeq_diabetes;
ODE_optionsG = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);
[~,xmodelG]   = ode15s(ODE_model, timespan, x0, ODE_optionsG,ptemp,c,timespan_i,datA_val,input);


E = [(xmodelG(:,2)-datGlu(:));
    (xmodelG(:,3)-datIns(:))*0.1];
