% comparison of free period and F-test for equispaced designs, choose one freq
% prior simulate signals of different freq and amplitude
clear
mode  = 'real';
n     = 12;
fmin  = 1;
fmax  = n/3;

addpath('../../MATLAB/utils')
switch mode
	case 'test'
		Nsamp = 5;
		Nperm = 10;
		nrep  = 10;
		Nfqf  = 10;
	case 'real'
		Nsamp = 5e2;
		Nperm = 1e3;
		nrep  = 3e3;
		Nfqf  = 5e2;
end


fprintf('Running n=%d\n',n)

tt =  (0:n-1)'/n;
if fmax==n/2
	cf=.95;
else
	cf=1;
end

Amps   = [1 2 3];
freqs  = 1:0.5:3;

fqf  = linspace(1,cf*fmax,Nfqf);
data = [];
tic;
for freq = freqs
	for Amp = Amps
		data_loc = NaN(nrep,5); 
		parfor rep=1:nrep
			acro            = 2*pi*rand;
			pwrFtest        = evalFtestPower(tt,freq,acro,Amp);
			pwrFree         = estFPPloc(tt,Nsamp,freq,acro,Amp,Nperm,fqf);
			data_loc(rep,:) = [Amp,freq,acro,pwrFtest,pwrFree];
		end
		data= [data;data_loc];
	end
end
toc
% write data to file
outFile = sprintf('mres_freePeriodSweep_n%d.csv', n);
writematrix(data, outFile);
fprintf('Saved results to %s\n', outFile);


function pwr = estFPPloc(tt,Nsamp,freq,acro,Amp,Nperm,fqf)
mu      = Amp*cos(2*pi*freq*tt -acro); % simulated signal
mu      = reshape(mu,[],1);
x       = mu + randn([length(tt),1,1,1,Nsamp]);

fqf     = reshape(fqf,1,1,[]);
[~,Q]   = getQuadForm(tt,fqf);
alpha   = .05;
pwrGrid = squeeze(fastMCTinfpower(Q,x,Nperm,alpha));
pwr     = min(pwrGrid,[],2);
end

