addpath('~/research/powerCHORD2/matlab/')
fmin    = 1;
fmax    = 16;
Nmeas   = 32;
Nfreq   = 1e3;

% differential evolution settings
settings.Npop       = 1e2;
settings.Niter      = 1e3;
settings.eps        = .01;
settings.useGPUglob = false;
settings.time_max    = Inf;


settings.method     = 'diffEvolve';
settings.diffEVsample = 'uniform';
[Tmat,eigfinal,scores] = pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);
fprintf('\n\nDiffEv FIT: %0.2f\n',max(eigfinal))

settings.method     = 'diffEvolveCR';
settings.CR         = .01;
[Tmat,eigfinal,scores] = pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);

fprintf('DiffEv CR:  %0.2f\n',max(eigfinal))

res = {Tmat,eigfinal,scores,settings};