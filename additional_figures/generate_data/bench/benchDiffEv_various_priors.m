% amplitude high or low (A=1 or A=2), try fmax=N/3 N/4 N/2 and show equispaced compared
% design obtained from optimizing the bound
addpath('../../../MATLAB/utils')

mode  = 'real';
Nmeas = 48;

Amps  = [1.5 2];
fmaxs = [Nmeas/4,Nmeas/3,Nmeas/2];

switch mode
    case 'tiny' % number of pvals = nrep*Nsamp 
		Nsamp  = 2;
		Nperm  = 2;
		Nfreq  = 2;
		Nacro  = 2;
		Nfq    = 2;
        nrep   = 2; 
	case 'test' 
		Nsamp  = 1e2;
		Nperm  = 1e2;
		Nfreq  = 24;
		Nacro  = 8;
		Nfq    = 48;
        nrep   = 1; 
	case 'real'
		Nsamp  = 50;
		Nperm  = 1e2;
		Nfreq  = 32;
		Nacro  = 16;
		Nfq    = 100;
        nrep   = 20; 
	case 'realhd'
		Nsamp  = 50;
		Nperm  = 5e2;
		Nfreq  = 32;
		Nacro  = 16;
		Nfq    = 500;
        nrep   = 20; 
end

fmin = 1;
DIFFEV = readtable('../../data/diffEvolveOutput.csv');
data_all=[];
for Amp=Amps
	fprintf('|On Amp %d\n',Amp)
    for fmax=fmaxs
		fprintf('||On fmax %d\n',fmax)
		
		% load diffEv design
		mt    = DIFFEV(DIFFEV.Nmeas==Nmeas & DIFFEV.fmin==fmin & DIFFEV.fmax == fmax,:);
		tt   = mt{:,9:(9+Nmeas-1)}';
		if length(tt)~=Nmeas
			error('check table')
		end    
        
        cf = 1;
		pwr      = estimateFreePeriodPower(tt,Nsamp,fmin,fmax,Nperm,cf,Amp,Nfreq,Nacro,Nfq,nrep);
        pwr      = reshape(pwr,[],1);
        data     = [Amp*ones(length(pwr),1) fmax*ones(length(pwr),1) linspace(fmin,fmax,Nfreq)' pwr];
        data_all = [data_all;data];
    end
end
outFile = sprintf('../../data/prF1c_diffEv_n%d_mode%s.csv', Nmeas,mode);
writematrix(data_all, outFile);
fprintf('Saved results to %s\n', outFile);
