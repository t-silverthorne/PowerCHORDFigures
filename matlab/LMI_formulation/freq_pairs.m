addpath('../utils/')

clear
pset=2;
switch pset
    case 1
        n     = 48; % clean pattern
        Nmeas = 12 ;
        freqs  = [1,10];
    case 2
        n     = 48; % genuine trade-off in lambda_min
        Nmeas = 12 ;
        freqs = [1,12];
end

[mu,eta]=run_yalmip(freqs,Nmeas,n);
tau   = ((1:n)/n - 1/n)';
plot(tau(value(mu)>1e-12),1,'.k')
mt = tau(value(mu)>1e-12);

[~,evals]=getMinEigMulti(mt,min(freqs),max(freqs),2,false,'all')
