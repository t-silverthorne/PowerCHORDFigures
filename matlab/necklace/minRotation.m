function r = minRotation(s)
n=length(s);
f=-1*ones(1,2*n);
k=0;
for j =1:2*n
    i = f(1+j-k-1);
    while i ~= -1 && s(1+mod(j,n))~= s(1+mod(k+i+1,n))
        if s(1+mod(j,n)) < s(1+mod(k+i+1,n))
            k = j-i-1;
        end
        i=f(1+i);
    end
    if i==-1 && s(1+mod(j,n)) ~= s(1+mod(k+i+1,n))
        if s(1+mod(j,n)) < s(1+mod(k+i+1,n))
            k=j;
        end
        f(1+j-k)=-1;
    else
        f(1+j-k)=i+1;
    end
end
r=s(1+mod(k:k+n-1,n));
end