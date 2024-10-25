%% 
nrep = 1e3;
%fprintf('Lemma 6.11: Tested %0.5g cases\n\n',1e9)
%% Lemma 6.9

tvec = NaN(1,nrep);
% (check series expression is correct)
for ii=1:nrep
    n1 = randsample(1:50,1);
    n2 = randsample(1:50,1);
    lambda = 10*rand(1);
    x = lambda*(0.5+rand(1));
    F_exact = ncfcdf(x,n1,n2,lambda);
    T = @(k,x,n1,n2,lambda) (lambda/2).^k*exp(-lambda/2).*...
        betainc(n1*x/(n2+n1*x),n1/2+k,n2/2)./factorial(k);
    F_approx = sum(T(0:50,x,n1,n2,lambda));
    tvec(ii) = abs(F_exact-F_approx)/F_exact; % check equal
end
max(tvec)
%%
% check beta identity is corect

tvec = NaN(1,nrep);
for ii=1:nrep
    b=1+10*rand(1); % integrals get difficult for a,b too small
    a=1+10*rand(1);
    z=rand(1);
    lhs=beta(a,b)*betainc(z,a,b) ;
    rhs=integral(@(t) t.^(a-1).*(1-t).^(b-1),0,z);
    tvec(ii)=abs(lhs-rhs)/rhs; % check equal
end
max(tvec)

%% Lemma 6.10
tvec1 = NaN(1,nrep);
tvec2 = NaN(1,nrep);
tvec3 = NaN(1,nrep);
for ii=1:nrep
    a=10*rand(1);
    b=10*rand(1);
    z=rand(1);
    % early identities
    tvec1(ii)=abs(beta(a,b)-gamma(a)*gamma(b)/gamma(a+b))/beta(a,b);
    tvec2(ii)=abs(betainc(z,a,b)- z^a*(1-z)^b/a/beta(a,b) - betainc(z,a+1,b))/betainc(z,a+1,b);
    lhs=betainc(z,a+1,b)/betainc(z,a,b);
    rhs=1 - (1-z)/(1 + (b-1)*integral(@(x) (x/z).^a.*(1-x).^(b-2),0,z )/(1-z)^(b-1));
    tvec3(ii)=abs(lhs-rhs)/lhs;
end
max(tvec1)
max(tvec2)
max(tvec3)

%% symbolic integral
syms x z b
assume(b>0)
assumeAlso(z>0)
assumeAlso(z<1)
int((1-x)^(b-2),x,0,z)

% limit
a=1e3
b=rand(1);
z=rand(1);
lhs=betainc(z,a+1,b)/betainc(z,a,b)
rhs=z

%% Lemma 6.11
tvec1 = NaN(1,nrep);
tvec2 = NaN(1,nrep);
tvec3 = NaN(1,nrep);
tvec4 = NaN(1,nrep);

for ii=1:nrep
    n1 = randsample(1:20,1);
    n2 = randsample(1:20,1);
    lambda = 10*rand(1);
    x = lambda*(0.5+rand(1));
    z = n1*x/(n2+n1*x);
    T = @(k,x,n1,n2,lambda) (lambda/2).^k.*exp(-lambda/2).*...
        betainc(n1*x/(n2+n1*x),n1/2+k,n2/2)./factorial(k);
    
    % check bound on series
    Iterm = @(k,x,n1,n2) betainc(n1*x/(n2+n1*x),n1/2+k,n2/2);
    k = randsample(1:50,1);
    tvec1(ii)=T(k+1,x,n1,n2,lambda)/T(k,x,n1,n2,lambda) <= lambda/2/(k+1);
    
    % check majorant for k>0
    prod(T(k,x,n1,n2,0:.1:100) <= Iterm(k,x,n1,n2)*k^k*exp(-k)/factorial(k))
    lhs=Iterm(k,x,n1,n2)*k^k*exp(-k)/factorial(k);
    rhs=T(k,x,n1,n2,2*k);
    tvec2(ii)=abs(lhs-rhs)/rhs;
    
    
    % finite difference derivative
    h = lambda*1e-3;
    (T(k,x,n1,n2,lambda+h)-T(k,x,n1,n2,lambda))/h;
    % get derivative symbolically
    lhs=abs(-0.5*T(k,x,n1,n2,lambda) + exp(-lambda/2)*lambda^(k-1)*...
        betainc(n1*x/(n2+n1*x),n1/2+k,n2/2)/factorial(k-1)/2^k);
    rhs=0.5*(T(k,x,n1,n2,lambda) + T(k-1,x,n1,n2,lambda));
    tvec3(ii) = lhs<rhs; % check derivative bound
    
    % majorant limit
    lhs = T(k+1,x,n1,n2,2*k+2)/T(k,x,n1,n2,2*k);
    rhs = Iterm(k+1,x,n1,n2)*exp(-1)*(1+1/k)^k/Iterm(k,x,n1,n2);
    tvec4(ii) = abs(lhs-rhs)/lhs;
end
prod(tvec1)
max(tvec2)
prod(tvec3)
max(tvec4)
%% Prop proof
tvec=NaN(1,nrep);
for ii=1:nrep
    N = randsample(5:50,1);
    n1 = randsample(1:20,1);
    n2 = randsample(1:20,1);
    lambda = 10*rand(1);
    h = lambda*1e-3;
    x = lambda*(0.5+rand(1));
    z = n1*x/(n2+n1*x);
    T = @(k,x,n1,n2,lambda) (lambda/2).^k.*exp(-lambda/2).*...
        betainc(n1*x/(n2+n1*x),n1/2+k,n2/2)./factorial(k);
    
    % syms lambda k
    % assume(k,'integer')
    % assumeAlso(lambda>0)
    % diff(exp(-lambda/2)*(lambda/2)^k/factorial(k),lambda)
    dT =@(k,x,n1,n2,lambda) ( ((1/2).^k.*k.*lambda.^(k - 1).*exp(-lambda/2))./factorial(k) - ((1/2).^k.*lambda.^k.*exp(-lambda/2))./(2*factorial(k))).*betainc(n1*x/(n2+n1*x),n1/2+k,n2/2);
    tvec(ii) = sum(dT(0:N,x,n1,n2,lambda)) <= -0.5*T(N,x,n1,n2,lambda);
end
prod(tvec)
%dT(N/2,x,n1,n2,lambda)
%(T(N/2,x,n1,n2,lambda+h) - T(N/2,x,n1,n2,lambda))/h
