%% Lemma 6.12
nrep = 1e7;
p = 3;%randsample(1:4,1);
q = 2;%randsample(1:p,1);
B = full(sprandsym(q,1));
B = B*B';
H = rand(q,p);

useH=true;
switch useH
    case true
        beta=rand(p,1,nrep);
        beta=pagemtimes(H,beta);
    case false
        beta=rand(q,1,nrep);
end        
        
beta = beta./pagenorm(beta);
est=min(pagemtimes(pagetranspose(beta),pagemldivide(B,beta)));
exact=1/max(eig(B));
exact<est
%% Lemma 6.13
nrep = 1e5;
n    = randsample(4:24,1);
mt   = rand(n,1,nrep);
freq = n*rand(1,1,nrep);
H    = [0 1 0; 0 0 1];
X    = [ones(n,1,nrep) cos(2*pi*freq.*mt) sin(2*pi*freq.*mt)];
tX   = [cos(2*pi*freq.*mt) sin(2*pi*freq.*mt)];
B    = pagemtimes(H,pagemldivide(pagemtimes(pagetranspose(X),X),H'));
bvec = pagetranspose(sum(tX,1));
lhs  = pageinv(B);
rhs  = pagemtimes(pagetranspose(tX),tX) - ...
    pagemtimes(bvec,pagetranspose(bvec))/n;

max(pagenorm(lhs-rhs))

%% Lemma 6.14

n = 5;
mt = reshape((1:n)/n -1/n,[],1);
X  = [ones(n,1) cos(2*pi*mt) sin(2*pi*mt)]
norm(X'*X - diag([n,n/2,n/2]))<1e-14
abs(sum(cos(2*pi*mt)))<1e-14
abs(sum(sin(2*pi*mt)))<1e-14



%% Lemma 6.15
syms theta
assume(theta,'real')
n = 5;
mt = reshape((1:n)/n -1/n,[],1);
tX  = [cos(2*pi*mt) sin(2*pi*mt)];
M  = tX'*tX/n;
E  = 0.5*eye(2);
f = [sin(theta);cos(theta)]
simplify(f'*E*f) <= min(eig(M))