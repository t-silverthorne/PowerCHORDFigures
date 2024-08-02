function genfix(t,p,n,d,a,b,Nmat)
global Nmat
if t >= d-1
    if mod(length(a)-1,p)==0 & max(a)==n-1 % only print necklaces
        uu = zeros(1,length(a)-1);
        uu(a(2:end))=1;
        Nmat(end+1,:)=uu;
    end
else
    tail = n-(d-t)+1;
    MAX  = a(t+2-p)+a(p+1);
    if MAX <= tail
        a(t+2)=MAX;
        b(t+2)=b(t+2-p);
        genfix(t+1,p,n,d,a,b);
        for i=(b(t+2)+1):1
            b(t+2) = i;
            genfix(t+1,t+1,n,d,a,b);
        end
        tail = MAX-1;
    end
    for j=tail:-1:(a(t+1)+1)
        a(t+2)=j;
        b(t+2)=1;
        genfix(t+1,t+1,n,d,a,b);
    end
end

