% amplitude high or low (A=1 or A=2), try fmax=N/3 N/4 N/2 and show equispaced compared
% design obtained from optimizing the bound
addpath('../../MATLAB/utils')

mode  = 'real';
Nmeas = 48;

Amps  = [1 2];
fmaxs = [Nmeas/4,Nmeas/3,Nmeas/2];
tu  = linspace(0,1,Nmeas+1);
tu  = tu(1:end-1)';

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
		Nsamp  = 1e2;
		Nperm  = 5e2;
		Nfreq  = 32;
		Nacro  = 32;
		Nfq    = 500;
        nrep   = 10; 
end

fmin = 1;
data_all=[];
for Amp=Amps
    for fmax=fmaxs
        if fmax==Nmeas/2
            cf = .95;
        else
            cf = 1;
        end
        pwr      = estimateFreePeriodPower(tu,Nsamp,fmin,fmax,Nperm,cf,Amp,Nfreq,Nacro,Nfq,nrep);
        pwr      = reshape(pwr,[],1);
        data     = [Amp*ones(length(pwr),1) fmax*ones(length(pwr),1) pwr];
        data_all = [data_all;data];
    end
end

outFile = sprintf('../data/prF1c_n%d_mode%s.csv', Nmeas,mode);
writematrix(data, outFile);
fprintf('Saved results to %s\n', outFile);
