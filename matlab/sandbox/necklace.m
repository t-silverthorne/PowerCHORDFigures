N = @(k,n) sum( arrayfun(@(ii) k^gcd(ii,n),1:n))/n
%%
n  = 72;
Nm = 12;
divs = divisors(gcd(Nm,n-Nm));

count=0;
for d=divs
    count = count + factorial(n/d)*eulerPhi(d)/factorial((n-Nm)/d)/factorial(Nm/d)/n;
end
count/10^4
%%
n  = 144;
Nm = 48;
A  = rand(n,n);
addpath('../')

t=rand(Nm,1e4);
tic
[~,f1]=getMinEigMulti(t,1,24,1e2,true,'min','fast');
toc
%%
tic
[~,f2]=getMinEigMulti(t,1,24,1e2,true,'min','slow');
toc
f1;
f2;


