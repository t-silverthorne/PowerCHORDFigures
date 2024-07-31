function [outputArg1,outputArg2] = genfix(t,p,n,d,a,b)
global Nmat
if t >= d-1
    Nmat(end+1,:)=a(2:n+1);
else
    tail = n-(d-t)+1;
    max  = a(t+2-p)+a(p+1);
    if max <= tail
        a(t+2)=max;
        b(t+2)=b(t+2-p);
        genfix(t+1,p,n,d,a,b);
        for i=(b(t+2)+1):1
            b(t+2) = i;
            genfix(t+1,t+1,n,d,a,b);
        end
        tail = max-1;
    end
    for j=tail:-1:(a(t+1)+1)
        a(t+2)=j;
        b(t+2)=i;
        genfix(t+1,t+1,n,d,a,b);
    end
end

