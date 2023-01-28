%%% example data (Leucine ingestion) used in the paper %%%
dataset = csvread('Example_data.csv');
inp = struct;
inp.OAAT.G.tim  = dataset(1:6,1)';
inp.OAAT.G.val  = dataset(1:6,2);
inp.OAAT.G.SE   = dataset(1:6,2).*0;
inp.OAAT.I.tim  = dataset(7:12,1)';
inp.OAAT.I.val  = dataset(7:12,2);
inp.OAAT.I.SE   = dataset(7:12,2).*0;
inp.OAAT.A.tim  = dataset(13:18,1)';
inp.OAAT.A.val  = dataset(13:18,2); % can be scaled to mmol/l if prefered %
inp.OAAT.A.SE   = dataset(13:18,2).*0; % can be scaled to mmol/l if prefered %

%%% Load parameter values, inputs and constants %%%
[p,c] = load_pHealthy_c([],inp.OAAT.G.val, inp.OAAT.I.val,1);
sim.Mg              = 0;   sim.Gpl             = 4.9;   sim.Ipl             = 6.7;    sim.Gint            = 0; sim.Id1=0;
data.Dmeal{1,1}     = 50e3;  data.t_meal_start   = 0;   data.Mb             = 70;            
k                   = 1;                   
timespan            = 0:0.1:150;     
[x0,input]          = load_x0_input(sim,data,k,timespan);

%%% Select parameters to estimate and the number of runs + random noise %%%
consideredP = [1 5 6 8 13]; runs = 25; settingI.rand = 25;

n_residuals = length(inp.OAAT.I.val) + length(inp.OAAT.G.val); n_ti= length(inp.OAAT.I.tim);    

%%% Interpolation Apl and dApl %%%
AA(:,1) = interp1([-20 -10 inp.OAAT.A.tim], [inp.OAAT.A.val(1) inp.OAAT.A.val(1) inp.OAAT.A.val'],timespan,'pchip');
AA(:,2) = gradient(AA(:,1)) ./ gradient(timespan(:));

for y = [5]
    z = y; 
    v = combnk(consideredP,z);
    residuals_multistart = zeros(1,runs,n_residuals,size(v,1));
    sim_glu = zeros(1,length(timespan),size(v,1));
    sim_ins = zeros(1,length(timespan),size(v,1));
    residuals = zeros(1,n_residuals,size(v,1));
    p_values = zeros(1,z,size(v,1));
    v_test = zeros(1,z,size(v,1));
    
    for comb=1:size(v,1)
        a=1;
        sim_glu_ = zeros(1,length(timespan));
        sim_ins_ = zeros(1,length(timespan));
        p_values_= zeros(1,z);
        residuals_= zeros(1,n_residuals);
        v_test_ = zeros(1,z);
        
        for i= 1 
        p0= p';
        namP= {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
        [p,c] = load_pHealthy_c([],inp.OAAT.G.val, inp.OAAT.I.val,1);
        [x0,input] = load_x0_input(sim,data,k,timespan); 
        x0(2)= inp.OAAT.G.val(1,1);  x0(3) = inp.OAAT.I.val(1,1);
        p(16) = inp.OAAT.G.val(i)';  p(17) = inp.OAAT.I.val(i)';  
        p0= p0(v(comb,:));
        sim.Gpl   = p(16);      sim.Ipl   = p(17);   
        namP= namP(v(comb,:));
        lbP= p0*0;  hbP= p0*1000;   
        pRest= p; options = optimset('Algorithm','trust-region-reflective','MaxFunEvals',500,'TolX',1e-30,'display','iter')

        p_init_rand= zeros(length(p0),runs);
        p_opt = p;
        Resvec = zeros(y,runs);
        Resvecnorm = zeros(1,runs);
        Resvecres = zeros(12,runs);
           
            for j = 1:runs             
                p_init_rand(:,j) = p0 + p0.*(rand(size(p0))/100*settingI.rand); 
               [Resvec(:,j) Resvecnorm(:,j) Resvecres(:,j)] = lsqnonlin(@errorFunction, p_init_rand(:,j), lbP, hbP, options,pRest,namP,inp.OAAT.G.tim,inp.OAAT.G.val,inp.OAAT.I.val,timespan,AA,x0,c,input,sim);
            end
            
            [M,I] = min(Resvecnorm);
            p_values_(a,:) = Resvec(:,I);
            residuals_(a,:) = Resvecres(:,I)';
            v_test_(a,:) = cell2mat(namP(1,:));
            
             for ind= 1:length(p0)
                    p_opt(namP{ind},1) = Resvec(ind,I);
             end 
            
            global t_saved Gpl_saved
            t_saved = 0;
            Gpl_saved = inp.OAAT.G.val(i,1)';
            ODE_model    = @diffeq_diabetes;
            ODE_optionsG = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);
            [sim.tGID,xmodelGID] = ode15s(ODE_model,timespan,x0,ODE_optionsG,p_opt,c,timespan,AA,input);

            sim_glu_(a,:) = xmodelGID(:,2)';
            sim_ins_(a,:) = xmodelGID(:,3)';        
            a=a+1;       
        end
        sim_glu(:,:,comb) = sim_glu_;
        sim_ins(:,:,comb) = sim_ins_;
        residuals(:,:,comb) = residuals_;
        p_values(:,:,comb) = p_values_;
        v_test(:,:,comb) = v_test_;
    end 
    R = struct()
    R.glu = sim_glu;
    R.ins = sim_ins;
    R.p = p_values;
    R.rs = residuals;
    R.multistart = residuals_multistart;
    R.vstart= v_test;
    filename = sprintf('PARAMS%02d.mat', y);
    save(filename, 'R') 
end

%% Results of the optimization %%
optis_5_Original = load('PARAMS05');

[p,c] = load_pHealthy_c([],inp.OAAT.G.val, inp.OAAT.I.val,1);
sim.Mg              = 0;    sim.Gpl             = p(16);   sim.Ipl             = p(17);    sim.Gint            = 0; 
data.Dmeal{1,1}     = 50e3;  data.t_meal_start   = 0;   data.Mb             = 70;      sim.Id1 = 0;       
k                   = 1;                   
timespan            = 0:0.1:150;     
[x0,input]          = load_x0_input(sim,data,k,timespan);

optis = optis_5_Original; 
optimal_params_values = optis.R.p;
optimal_params_number = optis.R.vstart;

namP = num2cell(optimal_params_number); 
p_0_opt = p;
for ind= 1:length(optimal_params_values)
       p_0_opt(namP{ind}) = optimal_params_values(ind);
end

global t_saved Gpl_saved
t_saved = 0;
Gpl_saved = sim.Gpl;
ODE_model    = @diffeq_diabetes;
ODE_optionsG = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

[sim.tGID,xmodelGID] = ode15s(ODE_model,timespan,x0,ODE_optionsG,p_0_opt,c,timespan,AA,input);

figure('Renderer', 'painters', 'Position', [10 10 400 800])
subplot(311)
plot(sim.tGID,AA(:,1)*0.001,'LineWidth',2.5,'Color', 'black'); %division by 1000 for umol/l --> mmol/l
hold on
e= errorbar(inp.OAAT.A.tim, inp.OAAT.A.val*0.001, inp.OAAT.A.SE*0.001,'x','Color','black','LineWidth',1.5);
ylabel('Plasma AA (mmol/L)')
ylim([1 4])
xlim([0 150])
H=gca;
H.LineWidth=2;
set(gca,'FontSize',13)

subplot(312)
plot(sim.tGID,xmodelGID(:,2),'LineWidth',2.5,'Color', 'red');
hold on
e= errorbar(inp.OAAT.G.tim, inp.OAAT.G.val, inp.OAAT.G.SE,'x','Color','black','LineWidth',1.5);
ylabel('Plasma Glucose (mmol/L)')
ylim([4 6])
xlim([0 150])
H=gca;
H.LineWidth=2;
set(gca,'FontSize',13)

subplot(313)
plot(sim.tGID,xmodelGID(:,3),'LineWidth',2.5,'Color', 'blue');
hold on
e= errorbar(inp.OAAT.I.tim, inp.OAAT.I.val, inp.OAAT.I.SE ,'x','Color','black','LineWidth', 1.5);
ylabel('Plasma Insulin (mU/L)')
ylim([0 75])
xlim([0 150])
H=gca;
H.LineWidth=2;
set(gca,'FontSize',13)

