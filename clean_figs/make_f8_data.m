nvec = [24,48,96]
cMat = NaN(3,17);

for ii=1:length(nvec)
    cMat(ii,:)=arrayfun(@(dd) enum(nvec(ii),dd),4:20)
end

writematrix(log10(cMat),'data/neck_counts.csv')

%%
plot(4:20,log10(cMat(1,:)))



function count=enum(n,d)
Nm = d;
divs = divisors(gcd(Nm,n-Nm));
count=0;
for d=divs
    count = count + factorial(n/d)*eulerPhi(d)/factorial((n-Nm)/d)/factorial(Nm/d)/n;
end
end
