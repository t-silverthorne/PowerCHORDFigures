Nm    = 48;
Npop  = 1e2;
fmin  = 1;
fmax  = 24;
Nfreq = 2^10;

t = rand(Nm,1);
useGPUglob=false;

Niter=100;
scores=NaN(Niter,1);
tic
for jj=1:Niter
    jj
    for ii=1:Nm
        tloc = t([1:(ii-1) (ii+1):Nm]);
        Jloc=@(s) -Out2(@getMinEigMulti,[s;tloc],fmin,fmax,Nfreq,useGPUglob);
        t(ii)=fminbnd(@(s) Jloc(s),0,1);
    end
    [~,scores(jj)]=getMinEigMulti(t,fmin,fmax,Nfreq,useGPUglob);
    
end
toc
close all
plot(scores,'.k')

function Z = Out2(FUN,varargin)
% Z = Out2(FUN,VARARGIN);
%
%	Provides the second output from the function
[~,Z] = FUN(varargin{:});
end