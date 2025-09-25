% compare wcp optimal vs equispaced for freq window [1,N/2] for various
% values of N, bandlimit the equispaced design to avoid numerical instability
clear;clf;
tic;
addpath('../../MATLAB/utils')
mode  = 'test';
Nmeas = 24;
cfact = .99;
Amp   = 2;
fmin  = 1;
fmax  = Nmeas/2;

% load WCP designs
data  = readtable('../data/diffEvolveOutput.csv');
mt    = data(data.Nmeas==Nmeas & data.fmin==fmin & data.fmax == fmax,:);

switch mode
	case 'test' % number of pvals = nrep*Nsamp 
		Nsamp  = 1e2;
		Nperm  = 1e2;
		Nfreq  = 24;
		Nacro  = 8;
		Nfq    = 48;
        nrep   = 1; 
	case 'real'
		Nsamp  = 1e2;
		Nperm  = 1e2;
		Nfreq  = 32;
		Nacro  = 32;
		Nfq    = 500;
        nrep   = 10; 
end

% define designs
nc  = 4; % number of clusters for comb design
tu  = linspace(0,1,Nmeas+1);
tu  = tu(1:end-1)';
tt_comb = make_comb(Nmeas,nc,.8/nc);
tt   = mt{:,9:(9+Nmeas-1)}';
if length(tt)~=Nmeas
	error('check table')
end    
tic;

% estimate free period power for equispaced design
fprintf('-----Running design 1 of 3-----\n')
pwr   = estimateFreePeriodPower(tt,Nsamp,fmin,fmax,Nperm,1,Amp,Nfreq,Nacro,Nfq,nrep);

fprintf('-----Running design 2 of 3-----\n')
pwru  = estimateFreePeriodPower(tu,Nsamp,fmin,fmax,Nperm,cfact,Amp,Nfreq,Nacro,Nfq,nrep);

fprintf('-----Running design 3 of 3-----\n')
pwrc  = estimateFreePeriodPower(tt_comb,Nsamp,fmin,fmax,Nperm,1,Amp,Nfreq,Nacro,Nfq,nrep);
freqs = linspace(fmin,fmax,Nfreq);

pwr  = reshape(pwr,[],1);
pwru = reshape(pwru,[],1);
pwrc = reshape(pwrc,[],1);

plot(pwr,'-k')
hold on
plot(pwru,'-b')
plot(pwrc,'-r')
toc
% outFile = sprintf('prF1c_n%d_Amp%d_mode%s.csv', Nmeas,Amp,mode);
% writematrix([pwru pwr pwrc], outFile);
% fprintf('Saved results to %s\n', outFile);
% toc
