function count=predictCount(n,d)
Nm = d;
divs = divisors(gcd(Nm,n-Nm));
count=0;
for d=divs
    count = count + factorial(n/d)*eulerPhi(d)/factorial((n-Nm)/d)/factorial(Nm/d)/n;
end

end

