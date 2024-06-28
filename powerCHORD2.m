% specify params
clear
freq  = 1.0;
Nmeas = 10;
A     = [0 0;
         1 0;
         0 1];
tau = 0:.01:1;
tau = tau(1:end-1);
nt  = length(tau);

c = reshape(cos(2*pi*freq*tau),[],1);
s = reshape(sin(2*pi*freq*tau),[],1);

dCS = [diag(c.*c) diag(c.*s);
       diag(c.*s) diag(s.*s) ];

CS=[c*c' s*c'
    c*s' s*s'];

IOmat = [ones(1,nt) zeros(1,nt);zeros(1,nt) ones(1,nt)];

Cmat = (dCS-CS/Nmeas);
Cm11 = Cmat(1:nt,1:nt);
Cm12 = Cmat(1:nt,(nt+1):2*nt);
Cm21 = Cmat((nt+1):2*nt,1:nt);
Cm22 = Cmat((nt+1):2*nt,(nt+1):2*nt);

mu =zeros(length(tau),1);
mu(randsample(1:length(tau),Nmeas,false)) = 1
[mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm21*mu mu'*Cm22*mu]


t = rand;
cvx_clear
%cvx_solver gurobi
cvx_begin
    variable mu(nt) binary;
    %minimize(mu'*(t*Cm11+ (1-t)*Cm22 + sqrt(t)*sqrt((1-t))*(Cm21+Cm12) )*mu)
    minimize(det(mu*mu'))
    subject to
        sum(mu)==Nmeas;
cvx_end



sanity_check =false;
if sanity_check
    mt = tau(mu>0);
    X     = [ones(length(mt),1) cos(2*pi*freq*mt)' sin(2*pi*freq*mt)'];
    B2 = A'*((X'*X)\A);
    inv(B2)
end