% Single frequency]
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
cvx_solver('mosek')
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
            % Mloc = [Nmeas cvec'*mu svec'*mu;
            %        [cvec'*mu;svec'*mu] Xr'*diag(mu)*Xr ];

            % M1=Xr'*diag(mu)*Xr;
            % v1=[cc 0 ; ss 0];
            % Mloc = diag([Nmeas,0,0]) + L1'*M1*L1 + L3'*v1*L2 + L2'*v1'*L3
            
            Mloc = Xr'*diag(mu)*Xr;
            lambda_min(Mloc) >= t;
        end
        sum(mu)==Nmeas;
        mu>=0;
cvx_end

%%
spec=[];spec_full=[];
for freq=freqs
    cc = cos(2*pi*freq*tau)'*mu;
    ss = sin(2*pi*freq*tau)'*mu;
    Xr = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
    Binv = Xr'*diag(mu)*Xr - [cc; ss]*[cc ss]/Nmeas;
    spec(end+1) = min(eig(Binv));
    spec_full(end+1) = min(eig([Nmeas cc ss;
             [cc;ss] Xr'*diag(mu)*Xr ]));
end
spec
spec_full


