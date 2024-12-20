addpath('../PowerCHORD/MATLAB/utils/')
addpath('../PowerCHORD/MATLAB/optim_methods/')

%% find first harmonic which cannot be resolved
Nvec = 6:48;
nh   = 100;
hvec = NaN(1,length(Nvec));
for ii=1:length(Nvec)
    N = Nvec(ii);
    mt = (1:N)/N -1/N;
    mt = mt';
    [freqs,eigs]=getMinEigMulti(mt,1,nh,nh,false,'all');
    bad_harm = find(abs((eigs-N/2)*2/N)>1e-8,1);
    if ~isempty(bad_harm)
        hvec(ii) = bad_harm;
    end
end
hvec
hvec

%% solve equivalent optimization problem and report power
ii             = 7
ovec = NaN(1,length(Nvec));
evec = NaN(1,length(Nvec));
for ii=1:length(Nvec)
    if Nvec(ii)>24
        [mu,eta,optim] = run_yalmip_timelim(1:hvec(ii),Nvec(ii),48*2,30);
    else
        [mu,eta,optim] = run_yalmip_timelim(1:hvec(ii),Nvec(ii),48,30);
    end
    ovec(ii)=optim.problem;
    evec(ii) = value(eta);
    ovec
    evec
end
%%
%%

function [mu,eta,optim]=run_yalmip_timelim(freqs,Nmeas,n,tlim)
tau   = ((1:n)/n - 1/n)';

eta   = sdpvar();
mu    = binvar(n,1);
F     = [sum(mu)==Nmeas,eta>=0,eta<=Nmeas/2];

for freq=freqs
    cvec = cos(2*pi*freq*tau);
    svec = sin(2*pi*freq*tau);
    uvec = ones(n,1);
    X=[uvec cvec svec];
    F=[F,X'*diag(mu)*X >= eta*eye(3)];
end
options=sdpsettings('solver','cutsdp');
options.cutsdp.maxtime=tlim;
optim=optimize(F,-eta,options)
end