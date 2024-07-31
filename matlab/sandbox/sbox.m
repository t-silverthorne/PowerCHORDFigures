
n  = 20;
Nm = 10;
mu = binvar(n,1);
nf = 10
sdpvar eta

F=sum(sum(mu))==Nm;

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

mu0 = [ones(Nm,1);zeros(n-Nm,1)];
eta0 = min(eig([mu0'*C11*mu0 mu0'*C12*mu0; mu0'*C12*mu0 mu0'*C22*mu0]));
warmstart(mu,mu0);
warmstart(eta,eta0);

optimize(F,-eta,sdpsettings('solver','bmibnb'))

%%
diag(A,A)