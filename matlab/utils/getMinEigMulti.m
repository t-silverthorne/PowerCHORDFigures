function [fv,em_fv] = getMinEigMulti(t,fmin,fmax,Nfreq,useGPU,returnType,method)
arguments % defaults
    t;
    fmin=1;
    fmax=24;
    Nfreq=1000;
    useGPU=true;
    returnType='min';
    method='fast';
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

switch method
    case 'fast'
        a = Cmat(1,1,:,:);
        b = Cmat(1,2,:,:);
        c = Cmat(2,1,:,:);
        d = Cmat(2,2,:,:);
        em_fv = a/2 + d/2 - (a.^2 - 2*a.*d + d.^2 + 4*b.*c).^(1/2)/2;
    otherwise
        if useGPU
            Cmat = gather(Cmat);
        end
        em_fv=min(pageeig(Cmat),[],1);
end
if strcmp(returnType,'min')
    em_fv = min(em_fv,[],3);
    em_fv = reshape(em_fv,1,[]);
end
end


