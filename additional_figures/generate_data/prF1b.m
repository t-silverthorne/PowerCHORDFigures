% for each harmonic, check how amplitude affects convergence of free period and
% F-test power, similar computation to fig1a, just different choice of which
% pars are randomised and which pars are stepped over discrete values

clear
addpath('../../MATLAB/utils')
Nsamp = 5e2;
Nperm = 1e3;
nrep  = 5e1;

n     = 12;
tt    = linspace(0,1,n+1);
tt    = tt(1:end-1)';

Namp  = 8;
Amin  = 1;
Amax  = 3;
Amps  = linspace(Amin,Amax,Namp);

fmin  = 1;
fmaxs = [n/2,n/3,n/4];
data_all = [];
for Amp=Amps
	for fmax=fmaxs
		if fmax==n/2
			cf=.95;
		else
			cf=1;
		end
		fqf  = linspace(fmin,cf*fmax,5e2);
		for freq=[1,n/4,n/3,n/2]
			data  = NaN(nrep,8);
			parfor rep=1:nrep
				acro        = 2*pi*rand;
				pwrFree     = estFPPloc(tt,Nsamp,freq,acro,Amp,Nperm,fqf);
				pwrFtest    = evalFtestPower(tt,freq,acro,Amp);
				vv          = [n fmin fmax Amp acro freq pwrFree pwrFtest];
				data(rep,:) = vv;
			end
			data_all = [data_all;data];
			writematrix(data_all, sprintf('data/results_prF1b_equi.csv'));
		end
	end
end

function [pwr,pwrGrid] = estFPPloc(tt,Nsamp,freq,acro,Amp,Nperm,fqf)
% wrapper for estimating power of free period model, calls fastMCTinfpower
% for actual power evaluation
mu      = Amp*cos(2*pi*freq*tt -acro); % simulated signal
mu      = reshape(mu,[],1);
x       = mu + randn([length(tt),1,1,1,Nsamp]);

fqf     = reshape(fqf,1,1,[]);
[~,Q]   = getQuadForm(tt,fqf);
alpha   = .05;
pwrGrid = squeeze(fastMCTinfpower(Q,x,Nperm,alpha));
pwr     = min(pwrGrid,[],2);
end

