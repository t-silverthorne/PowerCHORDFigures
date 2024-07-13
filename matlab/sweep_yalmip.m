mkdir 'yalmip_sweep'
fmin = 1:2:12;
fmax = 1:2:24;
df   = .5;
n    = 144;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,48]);

pars=[fmin(:) fmax(:) Nmeas(:)];


maxT=60;
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',maxT, ...
                    'bmibnb.maxiter',Inf);

ii      = str2double(getenv('SLURM_ARRAY_TASK_ID'));
fmin    = pars(ii,1);
fmax    = pars(ii,2);
Nmeas   = pars(ii,3);
fvec    = (fmin:df:fmax);


prob    = pCHORD(n,Nmeas,fvec,options);
res     = {fmin,fmax,Nmeas,prob};

save(fname,'res')