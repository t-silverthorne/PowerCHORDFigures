n       = 100;
N       = 12;
tau     = reshape((1:n)/n - 1/n,[],1);
mu      = reshape(randsample(1:n,N),[],1);
mub     = zeros(n,1);
mub(mu) = 1;
mt    = reshape(tau(mu),[],1);
freq  = 5*rand(1);

Xmt  = [ones(N,1) cos(2*pi*freq*mt) sin(2*pi*freq*mt)];
Xtau = [ones(n,1) cos(2*pi*freq*tau) sin(2*pi*freq*tau)];
norm(Xmt'*Xmt - Xtau'*diag(mub)*Xtau)

bvec = [cos(2*pi*freq*tau)'*mub;sin(2*pi*freq*tau)'*mub];
tX   = [cos(2*pi*freq*tau) sin(2*pi*freq*tau)]

Xschur = [N bvec';
          bvec tX'*diag(mub)*tX];

norm(Xschur-Xmt'*Xmt)


