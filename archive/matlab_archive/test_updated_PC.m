addpath('utils/')

fmin    = 1;
fmax    = 1;
Nmeas   = 6;
n       = 48;
fvec    = (fmin:1:fmax);

maxT=60*5;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf,...
                    'savesolveroutput',1,...
                    'warmstart',1, ...
                    'usex0',1);

method  = 'YALMIP';

settings.warmstart         = 'user'; % options: none, user, auto
settings.relax             = false;
settings.fix_gauge         = false;
settings.force_meas_region = true;
settings.no_meas_region    = true;
settings.FMRvec            = [zeros(n/2,1);1;zeros(n/4-1,1);1;zeros(n/4-1,1)];
settings.NMRvec            = [zeros(n/2,1);0;ones(n/4-1,1);0;ones(n/4-1,1)];
settings.FMRvec

settings.mu0    = [ones(Nmeas-sum(settings.FMRvec),1); zeros(sum(settings.FMRvec)+n-Nmeas,1)] +...
                            settings.FMRvec;

[prob,mu,eta,F] = pCHORD(n,Nmeas,fvec,options,method,settings);
mu
