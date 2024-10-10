%run('~/startup.m')


mkdir '~/research/powerCHORD2/matlab/sweep_yalmip_etabd'
addpath('~/research/powerCHORD2/matlab/')
addpath('~/research/powerCHORD2/matlab/utils/')
fmin = [1 2:2:12];
fmax = [1 2:2:24];
df   = .5;
n    = 144;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,40,48]);

pars=[fmin(:) fmax(:) Nmeas(:)];
maxT=60*6ww0;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf,...
                    'bmibnb.lpreduce',0,...
                    'savesolveroutput',1,...
                    'warmstart',1);

fmin    = 6;
fmax    = 20;
Nmeas   = 16;
fvec    = (fmin:df:fmax);

% fname = strcat('~/research/powerCHORD2/matlab/sweep_yalmip_etabd/', ...
%                'Nmeas_',num2str(Nmeas),...
%                '_fmin_',num2str(fmin),...
%                '_fmax_',num2str(fmax),...
%                '_maxT_',num2str(maxT),...
%                '_n_',num2str(n),...
%                '_lpred_','off',...
%                '.mat');
% 'fname defined'

[prob,mu,eta] = pCHORD_AC(n,Nmeas,fvec,options);

eta_v = value(eta)
mu_v  = value(mu)
res     = {fmin,fmax,Nmeas,prob,mu_v,eta_v};
save(fname,'res')

