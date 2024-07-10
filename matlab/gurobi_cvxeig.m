nm = 5;
nt = 40;
n=nt^2;
m=3;

a1=rand(n,m);
b1=rand(n,m);
c1=rand(n,m);

kcstr = 1e3;
Rcstr = randsample([0,1],kcstr*nt,true);
Rcstr = reshape(Rcstr,kcstr,[]);

R=[]
for ii=1:kcstr
    a = randsample(0:(nt-2),1);
    b = nt-2-a;
    R(ii,:) = [zeros(1,a*nt) Rcstr(ii,:) zeros(1,b*nt) -Rcstr(ii,:) ];
end
size(R)
cvx_setup

cvx_solver gurobi
cvx_begin
    
    variable s
    variable x(n) binary
    maximize(s)
    subject to 
        for ii=1:m
            s <= lambda_min([a1(:,ii)'*x b1(:,ii)'*x; b1(:,ii)'*x c1(:,ii)'*x])
        end
        % for ii=1:nt
        %     ind_start = (ii-1)*nt+1;
        %     ind_stop  = ii*nt;
        %     sum(x(ind_start:ind_stop )) <=nm
        % end
        R*x == 0
        sum(x)==nm*nm
        x>=0
        x<=1
cvx_end

rank(reshape(x,nt,nt))
nt