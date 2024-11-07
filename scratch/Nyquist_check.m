N  = 111;
mt = (1:N)/N - 1/N;
mt = mt';
freq=N/2;

X = [cos(2*pi*mt*freq) sin(2*pi*mt*freq)];

X'*X;

b = reshape(sum(X,1),2,[]);

X'*X-b*b'/N

X = [ones(N,1) cos(2*pi*mt*freq) sin(2*pi*mt*freq)];
H = [0 1 0; 0 0 1];

inv(H*((X'*X)\H'))
%%
min(eig(X'*X - b*b'/N))

%%
N=7
kvec=0:(N-1)
sum(exp(1j*pi*kvec)) + sum(exp(-1j*pi*kvec)) 