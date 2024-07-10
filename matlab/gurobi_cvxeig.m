nm = 5;
nt = 20;
n=nt^2;
m=1;

a1=rand(n,m);
b1=rand(n,m);
c1=rand(n,m);

cvx_setup
x0=rand(n,1);
cvx_solver gurobi
cvx_begin    
    variable s
    variable x(n) binary 
    maximize(s-norm(x-x0))
    subject to 
        for ii=1:m
            s <= lambda_min([a1(:,ii)'*x b1(:,ii)'*x; b1(:,ii)'*x c1(:,ii)'*x])
        end
        sum(x)==nm*nm
        X=[];
        for ii=1:nt
            istart = 1 + (ii-1)*nt;
            iend   = ii*nt;
            X=[X x(istart:iend)]
        end
            X==X';
        %ones(1,nt)*T*x ==nm
        x>=0;
        x<=1;
cvx_end

%M = nt*(nt-1)/2
%full(spdiags(ones(M,1),1,M,nt^2) -spdiags(ones(M,1),1+nt,M,nt^2))

rank(reshape(x,nt,nt))
nt