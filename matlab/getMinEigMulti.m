function [fv,em_fv] = getMinEigMulti(t,fmin,fmax,Nfreq,useGPU)
arguments % defaults
    t;
    fmin=1;
    fmax=24;
    Nfreq=1000;
    useGPU=true;
end

% unpack
Nm   = length(t);
fv   = linspace(fmin,fmax,Nfreq);
fv   = reshape(fv,1,1,[]);
%t    = reshape(t,[],1);

if useGPU
    t=gpuArray(t);
    fv=gpuArray(fv);
end

% build matrix
c    = cos(2*pi*fv.*t);
s    = sin(2*pi*fv.*t);

Cmat = [pagemtimes(pagetranspose(c),c) pagemtimes(pagetranspose(c),s);
      pagemtimes(pagetranspose(c),s) pagemtimes(pagetranspose(s),s)]- ...
      [sum(c,1).^2 sum(c,1).*sum(s,1); sum(c,1).*sum(s,1) sum(s,1).^2]/Nm;

if useGPU
    em_fv=min(pageeig(gather(Cmat)),[],1);
else
    em_fv=min(pageeig(Cmat),[],1);
end
end

