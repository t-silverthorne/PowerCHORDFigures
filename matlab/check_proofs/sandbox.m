p=3;
n=10;

X=rand(n,p);
eps=randn(n,1);
beta=rand(p,1);
y = X*beta + eps;
H = [0 1 0; 0 0 1];
B = H*((X'*X)\(H'));
Q = X*((X'*X)\(H'));

bhat = ((X'*X)\(X'*y));
num1 =  (H*bhat)'*(B\(H*bhat))
num2 =  y'*Q*((Q'*Q)\(Q'*y))
num3 =  (y-X*beta)'*Q*((Q'*Q)\(Q'*(y-X*beta)))

%%
Pu = X*((X'*X)\(X'));
den1 = y'*(eye(n)-Pu)*y
den2 = (y-X*beta)'*(eye(n)-Pu)*(y-X*beta)