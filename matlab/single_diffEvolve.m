addpath('~/research/powerCHORD2/matlab/')
fmin    = 1;
fmax    = 24;
Nmeas   = 48;
Nfreq   = 1e3;

% differential evolution settings
settings.method     = 'diffEvolve';
settings.Npop       = 1e3;
settings.Niter      = 1e4;
settings.eps        = .01;
settings.useGPUglob = true;

fname = strcat('diffev_Nmeas_',num2str(Nmeas),...
               '_fmin_',num2str(fmin),...
               '_fmax_',num2str(fmax),...
               '_fmax_',num2str(fmax),...
               '_Niter_',num2str(settings.Niter),...
               '_Npop_',num2str(settings.Npop),...
               '.mat');

res   = 1;
save(fname,'res');
[Tmat,eigfinal,scores] = pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);

res = {Tmat,eigfinal,scores,settings};
save(fname,'res');