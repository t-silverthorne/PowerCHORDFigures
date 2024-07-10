n=3;
uu=[]
for ii=1:1e2
    x=randsample([0,1],n^2,true);
    x=reshape(x,n,[]);
    %x=[1 0 1; 0 0 0;1 0 1];
    vals=[arrayfun(@(ii) sum(x(:,ii))==sum(x(ii,:)),1:n)];
    if all(vals) & sum(x)<n^2 & sum(x)>0
        uu(end+1)=rank(x);
        x
    end
end

%%
n=6
x=randsample([0,1],n,true);

x'*x