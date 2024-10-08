%run('~/startup.m')
addpath('utils/')
fmin    = 1;
fmax    = 12;
Nmeas   = 18;
n       = 48;
fvec    = linspace(fmin,fmax,2);

maxT=10;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf,...
                    'bmibnb.lpreduce',NaN,...
                    'savesolveroutput',1,...
                    'warmstart',1);

[prob,mu,eta,F] = pCHORD(n,Nmeas,fvec,options);
%%
check(F)
value(mu)
%%

fname = stretacat('~/research/powerCHORD2/matlab/yalmip_single/', ...
               'Nmeas_',num2str(Nmeas),...
               '_fmin_',num2str(fmin),...
               '_fmax_',num2str(fmax),...
               '_maxT_',num2str(maxT),...
               '_n_',num2str(n),...
               '.mat');
'fname defined'

'YALMIP completed'
eta_v = value(eta)
mu_v  = value(mu)
res     = {fmin,fmax,Nmeas,prob,mu_v,eta_v};
'res extracted'
save(fname,'res')
'file saved'