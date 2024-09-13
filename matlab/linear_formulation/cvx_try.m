
%% naive E-optimality
n= 100;
Nmeas=20;
tau = ((1:n)/n - 1/n)';

X = [ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)];

cvx_begin
    variables t 
    variable mu(n) binary
    maximize(t)
    subject to
        X'*diag(mu)*X>= t*eye(3,3)
        sum(mu)==Nmeas;
        mu>=0;
cvx_end

%%
n= 100;
Nmeas=20;
tau = ((1:n)/n - 1/n)';

X = [ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
H = [0 1 0; 0 0 1];
cvx_begin
    variable Z(3,2)
    variable mu(n) binary
    minimize(lambda_max(H*Z))
    subject to
        %H*Z             <= t*eye(2,2);
        sum(mu)         == Nmeas;
        X'*diag(mu)*X*Z == H';
        mu>=0;
cvx_end



%% power optimization
n     = 60;
Nmeas = 40;
tau   = ((1:n)/n - 1/n)';
fvec = [1,24];

cvx_begin
    variables b(2,length(fvec)) t
    variable mu(n) binary
    maximize(t)
    subject to
        for ff=1:length(fvec)
            freq=fvec(ff);
            Xr = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
            Xr'*diag(mu)*Xr - b(:,ff)*b(:,ff)'/Nmeas >= t*eye(2,2);
            cvec = cos(2*pi*freq*tau);
            svec = sin(2*pi*freq*tau);
            b(:,ff) == [cvec'*mu; svec'*mu];
        end
        sum(mu)==Nmeas;
        
        mu>=0;
        t>=0
        t<=Nmeas/2
        mu
cvx_end

plot(tau(mu>0),1,'.k')
