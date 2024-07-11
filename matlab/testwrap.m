
fvec=linspace(fmin,fmax,Nfreq)';
options = sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',60*10, ...
                    'bmibnb.maxiter',Inf);

pCHORD(n,Nm,fvec,options)