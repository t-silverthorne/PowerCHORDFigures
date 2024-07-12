clear all
addpath('../')

% roughly 10 minutes of compute time
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',60, ...
                    'bmibnb.maxiter',Inf, ...
                    'warmstart',1);


% harder problem: freq prior [1,24]
n       = 144;
Nm      = 48;
fvec    = linspace(1,24,2);
[prob,mu,eta]=pCHORD(n,Nm,fvec,options,true);

%% trial problem: single 24hr period with 8 measurements
n       = 24;
Nm      = 8;
options.bmibnb.maxtime=5;
freq    = rand*5;
fvec    = linspace(freq,freq,1);
[prob,mu,eta]=pCHORD(n,Nm,fvec,sdpsettings(),true);

%%
n       = 24;
Nm      = 8;
options.bmibnb.maxtime=5;
options.savesolveroutput=true;
options.savesolverinput=true;
freq    = 1;
fvec    = linspace(freq,freq,1);
[prob,mu,eta]=pCHORD(n,Nm,fvec,options);


%%
parfor ii=1:10
    print(ii)
end
