function [outputArg1,outputArg2] = gen(t,p,n,a)
global Nmat
if t>n
    Nmat(end+1,:)=a(2:n+1);
else
    a(t+1) = a(t+1-p);gen(t+1,p,n,a);
    for j=(a(t+1-p)+1):1
        a(t+1)=j;gen(t+1,t,n,a);
    end
end

end

