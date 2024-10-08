
clear
run('~/Documents/MATLAB/startup.m')


n=40;
Nm=10;
mu = binvar(n,1);

C11 = rand(n,n);
C22 = rand(n,n);
C12 = rand(n,n);
C11 = C11*C11';
C22 = C22*C22';
C12 = C12*C12';

F=sum(sum(mu))==Nm;

h=sum(mu);

options=sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpreduce',0,...
                    'bmibnb.maxtime',5);
A=[mu'*C11*mu mu'*C12*mu; mu'*C12*mu mu'*C22*mu];
mu0 = [ones(Nm,1);zeros(n-Nm,1)]
eta0 = min(eig([mu0'*C11*mu0 mu0'*C12*mu0; mu0'*C12*mu0 mu0'*C22*mu0]))
warmstart(mu,mu0)

sdpvar eta
warmstart(eta,eta0)
p2=optimize([sum(mu)==Nm,...
    sum(mu(1:n/2)) >= sum(mu((n/2+1):end)), ...
    sum(mu(1:n/4)) >= sum(mu((n/4+1):n/2)), ...
    A-eta*eye(2)>=0], ...
    -eta, ...
    sdpsettings('solver','bmibnb'))

% p1=optimize([sum(mu)==Nm], ...
%     -lambda_min(A), ...
%     sdpsettings('solver','moment'))
