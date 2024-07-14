mkdir 'sweep_diffevolve'

% differential evolution settings
settings.method     = 'diffEvolve';
settings.Npop       = 1e3;
settings.Niter      = Inf;
settings.time_max   = 10;
settings.eps        = .01;
settings.useGPUglob = true;

% define range to sweep over
fmin = 1:2:12;
fmax = 1:2:24;
Nfreq = 2^10;
n    = 144;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,48]);
pars=[fmin(:) fmax(:) Nmeas(:)];
pars=pars(pars(:,1)<=pars(:,2),:); % want fmin<=fmax

parpool(10)
parfor ii=1:size(pars,1)
	ii 
    fmin    = pars(ii,1);
    fmax    = pars(ii,2);
    Nmeas   = pars(ii,3);
		    
    fname = strcat('sweep_diffevolve/', ...
                   'Nmeas_',num2str(Nmeas),...
                   '_fmin_',num2str(fmin),...
                   '_fmax_',num2str(fmax),...
                   '_Niter_',num2str(settings.Niter),...
                   '_tmax_',num2str(settings.time_max),...
                   '.mat');
    [Tmat,eigfinal,scores] = pCHORDcts(Nmeas,fmin,fmax,Nfreq,settings);
    res=struct('Tmat',Tmat,...
               'eigfinal',eigfinal,...
               'scores',scores,...
               'settings',settings);
    save(fname,"-fromstruct",res)
end
