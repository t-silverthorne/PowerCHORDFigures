clear
n     = 24;
Nmeas = 8;
freq  = 1;

tau   = ((1:n)/n - 1/n)';

eta   = sdpvar();
mu    = binvar(n,1);
F     = [sum(mu)==Nmeas];

cvec = cos(2*pi*freq*tau);
svec = sin(2*pi*freq*tau);
cc   = cos(2*pi*freq*tau)'*mu;
ss   = sin(2*pi*freq*tau)'*mu;
Xr   = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
Mloc = [Nmeas cvec'*mu svec'*mu;
             [cvec'*mu;svec'*mu] Xr'*diag(mu)*Xr ];
F=[F,Mloc>=eta];

optimize(F,-eta)
plot(tau(value(mu)>0),1,'.k')
xlim([0,1])
value(eta)

%%
n     = 32;
Nmeas = 6;
freq  = 2;
tau   = ((1:n)/n - 1/n)';

eta   = sdpvar();
mu    = binvar(n,1);
F     = [sum(mu)==Nmeas];
cLHS  = 0;
for ii=1:n
    tt=tau(ii);
    cLHS = cLHS+[1 cos(2*pi*freq*tt) sin(2*pi*freq*tt);
                 cos(2*pi*freq*tt) cos(2*pi*freq*tt)^2 cos(2*pi*freq*tt)*sin(2*pi*freq*tt);
                 sin(2*pi*freq*tt) cos(2*pi*freq*tt)*sin(2*pi*freq*tt) sin(2*pi*freq*tt)^2]*mu(ii);
end
F=[F,cLHS>=eta]
options=sdpsettings('solver','sdpt3')
optimize(F,-eta,options)
plot(tau(value(mu)>0),1,'.k')
