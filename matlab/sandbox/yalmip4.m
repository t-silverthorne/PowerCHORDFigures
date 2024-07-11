clear all
addpath('../')

% roughly 10 minutes of compute time
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',60*10, ...
                    'bmibnb.maxiter',Inf);

%% trial problem: single 24hr period with 8 measurements
n       = 24;
Nm      = 8;
parfor ii=1:7
    freq = rand*5;
    fvec    = linspace(freq,freq,1);
    [mu,eta]=pCHORD(n,Nm,fvec,options);
end
%% harder problem: freq prior [1,24]
n       = 144;
Nm      = 48;
fvec    = linspace(1,24,24);
[mu,eta]=pCHORD(n,Nm,fvec,options);