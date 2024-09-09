'starting'
clear all
global Nmat
Nmat=[];
a(1)=0;

tic
gen(1,1,4,a);
toc

Nmat

%%
clear all
global Nmat
Nmat=[];
a=[0 0];

d = 6;
n = 144;
predictCount(n,d)
% expected number of fixed density binary necklaces

d=d+1;n=n+1;
jvals = ceil(n-d+1):-1:floor((n-1)/d+1);
tic
for j=jvals
    a(2) = j;
    b(2) = 1;
    genfix(1,1,n,d,a,b)
end
toc
size(Nmat)

fname = 'NeckMat_d_'+string(d);
save(fname,Nmat);


%%

function count=predictCount(n,d)
Nm = d;
divs = divisors(gcd(Nm,n-Nm));
count=0;
for d=divs
    count = count + factorial(n/d)*eulerPhi(d)/factorial((n-Nm)/d)/factorial(Nm/d)/n;
end
end
