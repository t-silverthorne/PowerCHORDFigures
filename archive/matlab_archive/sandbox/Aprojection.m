n = 1e2;
Nmeas = 12
tau = ((1:n)/n -1/n)';
freq = 5*rand;
X=[ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)]

cvx_solver Mosek
cvx_begin
    variable mu(n) binary
    maximize(lambda_min(A*A'*X'*diag(mu)*X*A*A'))
    subject to
        sum(mu)==Nmeas;
cvx_end
%%
n1=144;
n2=48;
tau = ((1:n1)/n1 -1/n1)';
freq1=1;
freq2=12;
A=[0 0;1 0;0 1];
X1=[ones(n1,1) cos(2*pi*freq1*tau) sin(2*pi*freq1*tau)];
X2=[ones(n1,1) cos(2*pi*freq2*tau) sin(2*pi*freq2*tau)];

s  = sdpvar(1);
mu = diag(binvar(n1,1));
F = [sum(sum(mu))==n2,s <= lambda_min(A*A'*X2'*mu*X2*A*A'), s<=lambda_min(A*A'*X1'*mu*X1*A*A')];
h = -s;
options=sdpsettings('solver','bnb');

optimize(F,h,options);

close all
plot(tau,diag(value(mu)),'.k')
%sol = value(mu)
%%
n1=24;
n2=8;
C11=rand(n1,n1);
C12=rand(n1,n1);
C22=rand(n1,n1);

mu = binvar(n1,1);

F=[sum(sum(mu))==n2];
h=-lambda_min([mu'*C11*mu mu'*C12*mu; mu'*C12*mu mu'*C22*mu]);
optimize(F,h)


%A=[0 0;1 0; 0 1];
