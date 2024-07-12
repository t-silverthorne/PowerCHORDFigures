run('~/startup.m')
mkdir '~/research/powerCHORD2/matlab/yalmip_sweep'
addpath('~/research/powerCHORD2/matlab/')
fmin = 1:2:12;
fmax = 1:2:24;
df   = .5;
n    = 144;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,48]);

pars=[fmin(:) fmax(:) Nmeas(:)];
pars=pars(pars(:,1)<=pars(:,2),:); % want fmin<=fmax

maxT=60*60;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf,...
                    'savesolveroutput',1);

ii      = str2double(getenv('SLURM_ARRAY_TASK_ID'));
fmin    = pars(ii,1);
fmax    = pars(ii,2);
Nmeas   = pars(ii,3);
fvec    = (fmin:df:fmax);

fname = strcat('~/research/powerCHORD2/matlab/yalmip_sweep/', ...
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