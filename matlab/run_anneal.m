
Nm=48;
fmin=1;
fmax=24;
Nfreq=1000;
cfun = @(t) -worstEig(t,fmin,fmax,Nfreq);
cfun(rand(48,1))

opts =optimoptions('simulannealbnd', ...
    'MaxIterations',1000, ...
    'HybridFcn','patternsearch',...
    'PlotFcn',{'saplotf','saplotbestf'})
t=simulannealbnd(cfun,rand(Nm,1),zeros(Nm,1),ones(Nm,1),opts);
close all
plot(t,1,'.k')


function Jt =worstEig(t,fmin,fmax,Nfreq)
[~,emins]=getMinEigMulti(t,fmin,fmax,Nfreq,false);
Jt=min(emins);
end