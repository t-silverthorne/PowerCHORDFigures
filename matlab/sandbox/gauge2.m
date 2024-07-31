clear
nf =1;
n=48;
m=12;
mu=binvar(n,1);
sdpvar eta

th1=(1:n)/n - 1/n;
obj = -eta

F = [sum(mu)==m];

for ii=1:nf
    C11 = rand(n,n);
    C22 = rand(n,n);
    C12 = rand(n,n);
    C11 = C11*C11';
    C22 = C22*C22';
    C12 = C12*C12';
    A=[mu'*C11*mu mu'*C12*mu; mu'*C12*mu mu'*C22*mu];
    
    F=[F,A-eta*eye(2)>=0];
end

options=sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpreduce',0,...
                    'bmibnb.maxtime',60);
prob1=optimize(F,obj,options);
v1 = value(obj);
mu0 = [ones(m/2,1); zeros(n-m,1); ones(m/2,1)];
% warmstart(mu,mu0);
% warmstart(eta,0);
F = [F,-2/n <= sum(diag(sin(2*pi*th1))*mu)/n <= 2/n];
prob2=optimize(F,obj,options);
v2 = value(obj);
v1
v2

prob1
prob2
%%


 -1/n <= sum(diag(sin(2*pi*th1))*mu)/n <= 1/n