addpath('../utils')
addpath('../optim_methods')

% optimize for 24hr and 2hr rhythms
freqs=[1 1.1]; % units of cycles/day
Nmeas=6;     % sample size
n    =24;     % grid coarsensss

% get optimal solution as binary vector
[mu,eta]=run_yalmip(freqs,Nmeas,n);
tau = (1:n)/n - 1/n;

% convert to time (units of days)
tau(value(mu)>0) % optimal solution

% convert to circadian time (units of hours)
24*tau(value(mu)>0)

% sanity check

tunif = (1:Nmeas)/Nmeas - 1/Nmeas
min(getMinEig(tau(value(mu)>0)',freqs(1)),getMinEig(tau(value(mu)>0)',freqs(2)))
min(getMinEig(tunif',freqs(1)),getMinEig(tunif',freqs(2)))