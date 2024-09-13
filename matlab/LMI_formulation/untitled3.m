close all
clear all
cvx_clear

L1 = [0 1 0;0 0 1];
L2 = [1 0 0;0 1 0];
L3 = [0 1 0;0 0 1];

n     = 32;
Nmeas = 6;
freqs = [1];
tau   = ((1:n)/n - 1/n)';
cvx_solver('gurobi')
cvx_begin
    variables t 
    variable mu(n) binary
    maximize(t)
    subject to
        for freq=freqs
            cvec = cos(2*pi*freq*tau);
            svec = sin(2*pi*freq*tau);
            cc = cos(2*pi*freq*tau)'*mu;
            ss = sin(2*pi*freq*tau)'*mu;
            Xr = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
            M  = [Nmeas cvec'*mu svec'*mu;
                   [cvec'*mu;svec'*mu] Xr'*diag(mu)*Xr ];
            M2 = [M(1,1) M(1,2);
                  M(2,1) M(2,2)];

            det_inv(M-t*eye(3)) >= 1;
            det_inv(M2-t*eye(2)) >= 1;
        end
        sum(mu)==Nmeas;
        mu>=0;
cvx_end
