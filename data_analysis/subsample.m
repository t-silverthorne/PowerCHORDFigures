addpath('../PowerCHORD/MATLAB/utils/')
addpath('../PowerCHORD/MATLAB/optim_methods/')

nvec = 8:48

% for n=nvec
%     if mod(48,n)==0
%         n
%     end
% end


Nmeas    = 14;
n        = 48;
freqs    = [1 2 3]*2;

tau = (1:n)/n - 1/n;

%% nmeas = 8
[mu,eta] = run_yalmip(freqs,8,n);
plot(tau(value(mu)>0)*48,1,'.k')
xlim([0,48])
eta
writematrix(value(mu),'mu_opt8.txt')
%% nmeas = 12
[mu,eta] = run_yalmip(freqs,12,n);
plot(tau(value(mu)>0),1,'.k')
writematrix(value(mu),'mu_opt12.txt')

%% nmeas = 16
[mu,eta] = run_yalmip(freqs,16,n);
plot(tau(value(mu)>0),1,'.k')
writematrix(value(mu),'mu_opt16.txt')
%% nmeas = 18
[mu,eta] = run_yalmip(freqs,18,n);
plot(tau(value(mu)>0),1,'.k')
writematrix(value(mu),'mu_opt18.txt')

%% nmeas = 24
[mu,eta] = run_yalmip(freqs,22,n);
plot(tau(value(mu)>0),1,'.k')

writematrix(value(mu),'mu_opt22.txt')
%%


%^value(eta)
Nmeas=48/8
mt = (1:Nmeas)/Nmeas - 1/Nmeas;

[getMinEig(mt',2),getMinEig(mt',4),getMinEig(mt',6)]
