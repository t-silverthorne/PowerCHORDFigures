function [fv,em_fv] = getMinEigMulti(t,fmin,fmax,Nfreq,useGPU,returnType)
arguments % defaults
    t;
    fmin=1;
    fmax=24;
    Nfreq=1000;
    useGPU=true;
    returnType='min';
end

% unpack
Nm   = size(t,1);
fv   = linspace(fmin,fmax,Nfreq);
fv   = reshape(fv,1,1,[]);

if size(t,2)>1
    t    = reshape(t,Nm,1,1,[]);
end

if useGPU
    t=gpuArray(t);
    fv=gpuArray(fv);
end

% build matrix
c    = cos(2*pi*fv.*t);
s    = sin(2*pi*fv.*t);

C1 = [pagemtimes(pagetranspose(c),c) pagemtimes(pagetranspose(c),s);
      pagemtimes(pagetranspose(c),s) pagemtimes(pagetranspose(s),s)];

C2 = [sum(c,1).^2 sum(c,1).*sum(s,1); sum(c,1).*sum(s,1) sum(s,1).^2];
Cmat =C1-C2/Nm;

if useGPU
    em_fv=min(pageeig(gather(Cmat)),[],1);
else
    em_fv=min(pageeig(Cmat),[],1);
end
if strcmp(returnType,'min')
    em_fv = min(em_fv,[],3);
    em_fv = reshape(em_fv,1,[]);
end
end

