function [mu,eta]=run_yalmip(freqs,Nmeas,n)
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
options=sdpsettings('solver','cutsdp');
optimize(F,-eta,options)
end