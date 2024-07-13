Nm    = 48;
fmin  = 1;
fmax  = 24;
Nfreq = 2^10;

settings.method     = 'diffEvolve';
settings.Npop       = 5e2;
settings.Niter      = 10;
settings.eps        = .01;
settings.useGPUglob = true;

[Tmat,eigfinal,scores] = pCHORDcts(Nm,fmin,fmax,Nfreq,settings)

