clear
n     = 48;
Nmeas = 20;
%freqs  = [1,12];
freqs  = [1,4];

tau   = ((1:n)/n - 1/n)';

eta   = sdpvar();
mu    = binvar(n,1);
F     = [sum(mu)==Nmeas,eta>=0,eta<=Nmeas/2];

for freq=freqs
    cvec = cos(2*pi*freq*tau);
    svec = sin(2*pi*freq*tau);
    uvec = ones(n,1);
    X=[uvec cvec svec];
    F=[F,X'*diag(mu)*X >= eta*eye(3)];
end
options=sdpsettings('solver','cutsdp')
optimize(F,-eta,options)
plot(tau(value(mu)>0),1,'.k')
xlim([0,1])
value(eta)