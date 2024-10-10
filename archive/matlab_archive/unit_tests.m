%% getMinEig
assert(abs(getMinEig([0 0.25 0.6 0.7]',0.3)-0.02552592)<1e-7);


freq=10*rand;
t = rand(25,1);
assert(abs(getMinEig(t,freq,'div')-getMinEig(t,freq,'inv'))<1e-12)

%% getMinEigMulti

% check agrees with getMinEig
t=rand(40,1);
[fv,em_fv]=getMinEigMulti(t,1,24,24);

for ii=1:length(fv)
    er(ii)=abs(getMinEig(t,fv(ii))-em_fv(ii));
end
assert(max(er)<1e-12)

% check using it on multiple measurement schedules agrees
Nm = 40;
t1=rand(Nm,1);
t2=rand(Nm,1);
[~,r1]=getMinEigMulti(t1,1,24,1000);
[~,r2]=getMinEigMulti(t2,1,24,1000);

tboth = [t1 t2];
tboth = reshape(tboth,Nm,1,1,[]);
[~,rboth] = getMinEigMulti(tboth,1,24,1000);

assert(all(rboth(:,:,:,1)==r1) & all(rboth(:,:,:,2)==r2))


%% getFourQuadBlocks
n = 24;
Nm= n-1;
freq = 5*rand;
tau= (1:n)/n -1/n;

t = tau(2:end);

mu = [0; ones(n-1,1)];

[Cm11,Cm12,Cm21,Cm22]=getFourQuadBlocks(n,Nm,freq);

min(eig([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm21*mu mu'*Cm22*mu]))

getMinEig(t',freq)
[~,fv]=getMinEigMulti(t',freq,freq,1,false)

