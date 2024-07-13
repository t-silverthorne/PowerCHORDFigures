run('~/startup.m')
mkdir '~/research/powerCHORD2/matlab/yalmip_single'
addpath('~/research/powerCHORD2/matlab/')
fmin    = 1;
fmax    = 24;
Nmeas   = 48;
n       = 144;
fvec    = (fmin:.5:fmax);

maxT=60*60*24;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf,...
                    'savesolveroutput',1,...
                    'warmstart',1);

fname = strcat('~/research/powerCHORD2/matlab/yalmip_single/', ...
               'Nmeas_',num2str(Nmeas),...
               '_fmin_',num2str(fmin),...
               '_fmax_',num2str(fmax),...
               '_maxT_',num2str(maxT),...
               '_n_',num2str(n),...
               '.mat');
'fname defined'

[prob,mu,eta] = pCHORD(n,Nmeas,fvec,options);
'YALMIP completed'
eta_v = value(eta)
mu_v  = value(mu)
res     = {fmin,fmax,Nmeas,prob,mu_v,eta_v};
'res extracted'
save(fname,'res')
'file saved'