clear all
addpath('../')
n    = 144;
Nm   = 48;
tau  = ((1:n)/n -1/n)';
fvec = linspace(1,24,24);
sdpvar eta;
mu = binvar(n,1);

F=sum(sum(mu))==Nm;

for freq=fvec
    [Cm11,Cm12,~,Cm22] = getFourQuadBlocks(n,Nm,freq);
    F = [F, eta <=lambda_min([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm12*mu mu'*Cm22*mu])];
end
h =-eta;

options=sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','mosek', ... % can be gurobi or mosek
                    'bmibnb.maxtime',10);
%warning('off','all')
optimize(F,h,options)
%warning('on','all')