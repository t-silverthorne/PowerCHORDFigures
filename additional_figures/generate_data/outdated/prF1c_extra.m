% compare wcp optimal vs equispaced for freq window [1,N/2] for various
% values of N, bandlimit the equispaced design to avoid numerical instability
clear;clf;
tic;
mode  = 'real';
Nmeas = 48;
cfact = .95;
Amp   = 5;
fmin  = 1;
fmax  = Nmeas/2;

% load designs
data  = readtable('diffEvolveOutput.csv');
mt    = data(data.Nmeas==Nmeas & data.fmin==fmin & data.fmax == fmax,:);

addpath('../../MATLAB/utils')
switch mode
    case 'nothing'
        Nsamp  = 2;
		Nperm  = 2;
		Nfreq  = 2;
		Nacro  = 3;
		Nfq    = 5;
        nrep   = 3; 
	case 'test' % number of pvals = nrep*Nsamp 
		Nsamp  = 1e2;
		Nperm  = 1e2;
		Nfreq  = 24;
		Nacro  = 32;
		Nfq    = 500;
        nrep   = 10; 
	case 'real'
		Nsamp  = 1e2;
		Nperm  = 1e2;
		Nfreq  = 32;
		Nacro  = 16;
		Nfq    = 500;
        nrep   = 10; 
	case 'real-hd'
		Nsamp  = 1e2;
		Nperm  = 1e3;
		Nfreq  = 32;
		Nacro  = 16;
		Nfq    = 500;
        nrep   = 100; 
end

% estimate free period power for permutation optimised design
fig  = openfig('optim_results/test_Nmeas48_Amp10_MaxIter100_fmin1_fmax24_Nfreqch64_Nacroch64_Nsampch100_NfqTinf64_NfqT25000_Nsampmc500_Nfreqmc64_Nacromc64_Npermmc1000.fig', ...
               'invisible');
ax   = findall(fig,'type','axes');
ax_sorted = flipud(ax);
top_ax = ax_sorted(3);        % top panel
dots = findobj(top_ax,'Type','Line','Marker','.');
x    = get(dots,'XData');
y    = get(dots,'YData');
x    = x(:);
tt   = cell2mat(x(:));
tt   = reshape(tt,[],1);

if length(tt)~=Nmeas
	error('check table')
end    
tic;

% estimate free period power for equispaced design
pwr     = estimateFreePeriodPower(tt,Nsamp,fmin,fmax,Nperm,1,Amp,Nfreq,Nacro,Nfq,nrep);
freqs   = linspace(fmin,fmax,Nfreq);
pwr     = reshape(pwr,[],1);
outFile = sprintf('prF1c_extra_n%d_Amp%d_mode%s.csv', Nmeas,Amp,mode);
writematrix(pwr, outFile);
fprintf('Saved results to %s\n', outFile);
toc
% plot(squeeze(freqs),squeeze(pwr),'-k')
% hold on
% plot(squeeze(freqs),squeeze(pwru),'-b')
% drawnow
% ylim([0,1])
% toc
