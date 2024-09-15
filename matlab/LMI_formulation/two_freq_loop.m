addpath('../utils/')

n=48;
Nmeas=12;

freqs=[2 4 6 8 10 12];

Smat = NaN(length(freqs),n);
tic
for ii=1:length(freqs)
    [mu,eta]=run_yalmip([1 freqs(ii)],Nmeas,n);
    Smat(ii,:)=value(mu);
end
X = zeros(size(Smat));
X(Smat>1e-12)=1;
X
writematrix(X,'cutsdp_sols.csv')
%%