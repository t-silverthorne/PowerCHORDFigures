addpath('../utils/')
addpath('../optim_methods/')
testing = false;
dirname = 'diffEvolveOutput/';
mkdir(dirname); 

% differential evolution settings
settings.method     = 'diffEvolveCR';
settings.CR         = .05;
settings.Npop       = 1e3;
settings.Niter      = Inf;
settings.eps        = .05;
settings.useGPUglob = true;
switch testing
    case true
        settings.time_max   = 10;
    otherwise
        settings.time_max   = 60*60;
end


% define range to sweep over
fmin = [1 2:2:12];
fmax = [1 2:2:24];
Nfreq = 2^10;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,40,48]);
pars=[fmin(:) fmax(:) Nmeas(:)];
pars=pars(pars(:,1)<=pars(:,2),:); % want fmin<=fmax

switch testing
    case true
        ii = 200;
    otherwise
        ii=str2double(getenv('SLURM_ARRAY_TASK_ID'));
end


fmin    = pars(ii,1);
fmax    = pars(ii,2);
Nmeas   = pars(ii,3);
        
fname = strcat(dirname,...
                'Nmeas_',num2str(Nmeas),...
                '_fmin_',num2str(fmin),...
                '_fmax_',num2str(fmax),...
                '_Niter_',num2str(settings.Niter),...
                '_tmax_',num2str(settings.time_max),...
                '.mat');
[Tmat,eigfinal,scores] = diffEvolve(Nmeas,fmin,fmax,Nfreq,settings);
res=struct('Tmat',Tmat,...
            'eigfinal',eigfinal,...
            'scores',scores,...
            'settings',settings);
save(fname,"-fromstruct",res)
[turner97@beluga1 cluster_compute_gpu]$ :q
[mii] :q not found! Similar commands: "sq", "sqm", "mb"
[turner97@beluga1 cluster_compute_gpu]$ cat runDiffEvolve.m 
addpath('../utils/')
addpath('../optim_methods/')
testing = false;
dirname = 'diffEvolveOutput/';
mkdir(dirname); 

% differential evolution settings
settings.method     = 'diffEvolveCR';
settings.CR         = .05;
settings.Npop       = 1e3;
settings.Niter      = Inf;
settings.eps        = .05;
settings.useGPUglob = true;
switch testing
    case true
        settings.time_max   = 10;
    otherwise
        settings.time_max   = 60*60;
end


% define range to sweep over
fmin = [1 2:2:12];
fmax = [1 2:2:24];
Nfreq = 2^10;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,24,32,40,48]);
pars=[fmin(:) fmax(:) Nmeas(:)];
pars=pars(pars(:,1)<=pars(:,2),:); % want fmin<=fmax

switch testing
    case true
        ii = 200;
    otherwise
        ii=str2double(getenv('SLURM_ARRAY_TASK_ID'));
end


fmin    = pars(ii,1);
fmax    = pars(ii,2);
Nmeas   = pars(ii,3);
        
fname = strcat(dirname,...
                'Nmeas_',num2str(Nmeas),...
                '_fmin_',num2str(fmin),...
                '_fmax_',num2str(fmax),...
                '_Niter_',num2str(settings.Niter),...
                '_tmax_',num2str(settings.time_max),...
                '.mat');
[Tmat,eigfinal,scores] = diffEvolve(Nmeas,fmin,fmax,Nfreq,settings);
res=struct('Tmat',Tmat,...
            'eigfinal',eigfinal,...
            'scores',scores,...
            'settings',settings);
save(fname,"-fromstruct",res)

