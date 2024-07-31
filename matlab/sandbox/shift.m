n=500
m=100
inds=randsample(n,m)
mu=zeros(n,1)
mu(inds)=1
mu
w=(1:n)'
%%
score=NaN(1,n)

for ii=1:n
    mu = [mu(end); mu(1:end-1)]
    score(ii) = w'*mu
end
plot(w,score,'.k')

%%
score=NaN(1,n)
th1=(1:n)/n - 1/n;

for ii=1:n
    mu = [mu(end); mu(1:end-1)];
    score(ii) =  sum(diag(sin(2*pi*th1))*mu)/n ;
end
plot(w,score,'.k')