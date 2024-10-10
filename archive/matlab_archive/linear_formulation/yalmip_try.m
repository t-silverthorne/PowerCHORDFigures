clear
%% Single frequency YALMIP using BMI solver

n     = 24;
Nmeas = 8;
freq  = 1;
tau   = ((1:n)/n - 1/n)';
Z     = sdpvar(3,2);
eta   = sdpvar();
mu    = binvar(n,1);
X     = [ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
H     = [0 1 0; 0 0 1];

F     = [X'*diag(mu)*X*Z==H', H*Z <= eta,sum(mu)==Nmeas,eta>=2/Nmeas];
optimize(F,eta)

%% Multi-frequency YALMIP using Gurobi

n     = 72;
Nmeas = 40;
fvec  = [1,12,24];
tau   = ((1:n)/n - 1/n)';
Z     = sdpvar(3,2,length(fvec));
eta   = sdpvar();
mu    = binvar(n,1);

F     = [sum(mu)==Nmeas,eta>=2/Nmeas];
H     = [0 1 0; 0 0 1];
for ii=1:length(fvec)
    freq  = fvec(ii);
    Zloc  = Z(:,:,ii);
    X     = [ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
    F     = [F,X'*diag(mu)*X*Zloc==H', H*Zloc <= eta];
end
options=sdpsettings('solver','gurobi')
optimize(F,eta,options)

%%
%%
Xr = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)];

F  = [Xr'*diag(mu)*Xr - ]
%%


K = sdpvar(1,2);
P = sdpvar(2);
A = [-1 2;2 -3];
B = [1;1];
MyLMI = [ [P (A + B*K)';(A + B*K) inv(P)] >= 0];
optimize(MyLMI)
