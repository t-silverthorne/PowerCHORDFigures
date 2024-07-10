Nm=8;
freqs='single'
switch freqs
    case 'multi'
        fmin=1;
        fmax=24;
        Nfreq=1000;
    case 'single'
        fmin=1;fmax=1;Nfreq=1;
end

cfun = @(t) -worstEig(t,fmin,fmax,Nfreq);
cfun(rand(48,1))

opts =optimoptions('simulannealbnd', ...
    'MaxIterations',1000, ...
    'HybridFcn','patternsearch',...
    'PlotFcn',{'saplotf','saplotbestf'})
tic;
t=simulannealbnd(cfun,rand(Nm,1),zeros(Nm,1),ones(Nm,1),opts);
toc
close all
plot(t,1,'.k')


function Jt =worstEig(t,fmin,fmax,Nfreq)
[~,emins]=getMinEigMulti(t,fmin,fmax,Nfreq,false);
Jt=min(emins);
end