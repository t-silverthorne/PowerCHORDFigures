cd('~/research/powerCHORD2/matlab')
fmin  = 1;
fmax  = 24;
Nfreq = 1e3;
Nmeas = 32 

% try diffEvolve
settings.useGPUglob = true;
settings.method     = 'diffEvolve';
settings.Niter      = 1e3;
settings.eps        = .01;
settings.Npop       = 500;

tic;
[u,ef,~]=pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);
max(ef)
Tmax=toc 

% try simulated annealing 
settings.useGPUglob = false;
settings.method     = 'simulAnneal';
settings.Niter      = 100;
settings.maxTime    = Tmax;

tic;
[u,ef,~]=pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);
ef
toc
