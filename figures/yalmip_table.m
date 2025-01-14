addpath('../PowerCHORD/MATLAB/utils/')
addpath('../PowerCHORD/MATLAB/optim_methods/')

% optimize for 24hr and 2hr rhythms
freqs=[1 12]; % units of cycles/day
Nmeas=12;     % sample size
n    =48;     % grid coarsensss

% get optimal solution as binary vector
[mu,eta]=run_yalmip(freqs,Nmeas,n);
tau = (1:n)/n - 1/n;

% convert to time (units of days)
tau(value(mu)>0) % optimal solution

% convert to circadian time (units of hours)
24*tau(value(mu)>0)

% sanity check
value(eta)==6