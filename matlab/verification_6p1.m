%% Lemma 6.1 (rank condition in terms of distinct phase)
N   = randsample(4:10,1);
phi = [2*pi*rand(3,1);zeros(N-3,1)];
length(phi)==N;
X   = [ones(N,1) cos(phi) sin(phi)];
rank(X)==3

%% Lemma 6.2 (positivity)
N   = randsample(4:20,1);
p   = randsample(1:N,1);
X = rand(N,p);
all(eig(X'*X)>0)

%% Lemma 6.3 (positivity with random hypothesis matrix)
N   = randsample(4:20,1);
p   = randsample(1:N,1);
q   = randsample(1:p,1);
A   = full(sprandsym(p,1,rand(1,p)))
H   = rand(q,p);
all(eig(H*A*H')>0)

%% Lemma 6.4 (check nullspace of projector)
N    = randsample(4:15,1);
p    = randsample(1:N,1);
q    = randsample(1:p,1);
X    = rand(N,p);
H    = rand(q,p);

Q    = X*((X'*X)\(H'));
Pell = Q*((Q'*Q)\(Q'));
size(Pell)
nb   = null(H)
max(abs(Pell*X*nb))

%% Cor 6.4.1 - specific to cosinor model (check equality of projectors)
N    = randsample(4:15,1);
p    = 3;%randsample(1:N,1);
q    = 2;randsample(1:p,1);
X    = [ones(N,1) rand(N,p-1)];
H    = [0 1 0;0 0 1];

Pr =  ones(N,N)/N
Pu =  X*((X'*X)\X');
Q  =  X*((X'*X)\(H'));
Pl =  Q*((Q'*Q)\(Q'));
max(max(abs(Pl - (Pu-Pr))))


%% Lemma 6.6 (null distribution of F statistic)
nrep=1e4;
N  = randsample(10:25,1);
p  = randsample(3:4,1);
q  = randsample(1:(p-1),1);
X  = rand(N,p);
H  = rand(q,p);

nb    = null(H);
cvec  = rand(size(nb,2),1);
beta  = nb*cvec;

sigma = 3*rand(1);
Ydat  = X*beta + sigma*randn(N,1,nrep);

Pu       = X*((X'*X)\X');
bhat     = pagemtimes((X'*X)\X',Ydat);
Hbhat    = pagemtimes(H,bhat);

num   = reshape( pagemtimes(pagetranspose(Hbhat),pagemldivide( H*((X'*X)\H'),Hbhat)),[nrep,1,1]);
denom = reshape(pagemtimes(pagetranspose(Ydat),pagemtimes(eye(N)-Pu,Ydat)),[nrep,1,1]);
Fdat  = (N-p)*num./denom/(q);
close all
histogram(Fdat,floor(sqrt(nrep)),'Normalization','pdf')
hold on
x = 0:.01:max(Fdat);
plot(x,fpdf(x,q,N-p))

%% Lemma 6.7 (independence of disjoint projectors)
nrep = 1e4;
n    = 2*randsample(10:20,1);
X    = randn(n,1,nrep);

A = full(sprandsym(n,1));

[V,~]=eig(A);

Q1 = V(:,(1:(n/2)));
Q2 = V(:,((n/2+1):n));

P1 = Q1*((Q1'*Q1)\Q1');
P2 = Q2*((Q2'*Q2)\Q2');

max(max(abs(P1^2-P1)));
max(max(abs(P2^2-P2)));

max(max(abs(P1*P2)))

stat1 = reshape(pagenorm(pagemtimes(P1,X),2),nrep,1,1);
stat2 = reshape(pagenorm(pagemtimes(P2,X),2),nrep,1,1);

clf
plot(stat1,stat2,'.k')
cov(stat1,stat2)
%% Lemma 6.8
nrep = 1e4;
N  = randsample(10:25,1);
p  = randsample(3:6,1);
q  = randsample(1:(p-1),1);
X  = rand(N,p);
H  = rand(q,p);

nb    = null(H);
cvec  = rand(size(nb,2),1);
beta  = H'*rand(q,1);

mu    = 0*rand(N,1);
sigma = 3*rand(1);
Ydat  = X*beta + mu+ sigma*randn(N,1,nrep);

B   = H*((X'*X)\H');
ncp = (H*beta)'*(B\(H*beta))/sigma^2;

Pu       = X*((X'*X)\X');
bhat     = pagemtimes((X'*X)\X',Ydat);
Hbhat    = pagemtimes(H,bhat);

num   = reshape( pagemtimes(pagetranspose(Hbhat),pagemldivide( H*((X'*X)\H'),Hbhat)),[nrep,1,1]);
denom = reshape(pagemtimes(pagetranspose(Ydat),pagemtimes(eye(N)-Pu,Ydat)),[nrep,1,1]);
Fdat  = (N-p)*num./denom/(q);
close all
histogram(Fdat,floor(sqrt(nrep)),'Normalization','pdf')
hold on
x = 0:.01:max(Fdat);
plot(x,ncfpdf(x,q,N-p,ncp))