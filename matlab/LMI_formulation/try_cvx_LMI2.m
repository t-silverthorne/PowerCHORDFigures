%% Check matrices agree
clear

n    = 100;
N    = 11;
mu   = [ones(N,1);zeros(n-N,1)];
tau  = ((1:n)/n - 1/n)';
freq = 1;
cvec = cos(2*pi*freq*tau);
svec = sin(2*pi*freq*tau);
Xr   = [cvec svec]
ovec = zeros(n,1);
L1 =[ovec cvec svec]'

[0 0; 1 0; 0 1]*Xr'*diag(mu)*Xr*[0 0; 1 0; 0 1]'

N*diag([1,0,0])+ ...
L1*diag(mu)*ones(n,3)*diag([1,0,0]) + ...
    (L1*diag(mu)*ones(n,3)*diag([1,0,0]))' + ...
    [0 0; 1 0; 0 1]*Xr'*diag(mu)*Xr*[0 0; 1 0; 0 1]'


[N cvec'*mu svec'*mu;
 cvec'*mu cvec.^2'*mu (cvec.*svec)'*mu;
 svec'*mu (cvec.*svec)'*mu svec.^2'*mu ]

uvec = ones(n,1);
[uvec cvec svec]'*diag(mu)*[uvec cvec svec]
%%
m=3
addpath('.')
cvx_setup
n    = 500;
N    = 40;
mu   = [ones(N,1);zeros(n-N,1)];
tau  = ((1:n)/n - 1/n)';
freqs = 1:.1:24
cvec = cos(2*pi*freq*tau);
svec = sin(2*pi*freq*tau);
uvec = ones(n,1);

cvx_solver('mosek')
m=3;
cvx_begin sdp
    variables t 
    variable mu(n)  
    maximize(t)
    subject to
        uvec = ones(n,1);
        for freq=freqs
            cvec = cos(2*pi*freq*tau);
            svec = sin(2*pi*freq*tau);
            X=[uvec cvec svec];
            X'*diag(mu)*X >= t*eye(m);
        end
        sum(mu)==N;
        mu>=0;
        mu<=1;
cvx_end

plot(tau,mu)