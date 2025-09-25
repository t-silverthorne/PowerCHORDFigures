% compare Tinf, T2, T2bnd, and exact power estimation
clear;
addpath('../../MATLAB/utils/');

nvals  = 100:200:1100;
dat    = [];
tmeth  = 'tictoc';
n = 24;
for Nperm=nvals
    fprintf('Running Nperm = %d\n',Nperm)
    freq  = rand;
    acro  = 2*pi*rand;
    Amp   = rand;
    tt    = rand(n,1);

    ntrep = 100; % number of times to benchmark, currently unused

    Nsamp = 1e3;
    fqf   = reshape(linspace(1,2,1e3),1,1,[]);
    [Q2,Qinf] = getQuadForm(tt,fqf);
    wrap_Tinf = @(tt,freq,acro,Amp) estTinfloc(tt,Nsamp,freq,acro,Amp,Nperm,Qinf);
    wrap_T2   = @(tt,freq,acro,Amp) estT2loc(tt,Nsamp,freq,acro,Amp,Nperm,Q2);
    wrap_T2b  = @(tt,freq,acro,Amp) Jfun(tt,freq,acro,Q2,Amp,Nsamp);

    fprintf('  timing 1\n')
    [mm1,iqr1] = repTimeFunc(@evalFtestPower,ntrep,tmeth,tt,freq,acro,Amp);

    fprintf('  timing 2\n')
    [mm2,iqr2] = repTimeFunc(wrap_Tinf,ntrep,tmeth,tt,freq,acro,Amp);

    fprintf('  timing 3\n')
    [mm3,iqr3] = repTimeFunc(wrap_T2,ntrep,tmeth,tt,freq,acro,Amp);

    fprintf('  timing 4\n')
    [mm4,iqr4] = repTimeFunc(wrap_T2b,ntrep,tmeth,tt,freq,acro,Amp);
    datloc =   [Nperm mm1 iqr1(:)' 1;
                Nperm mm2 iqr2(:)' 2;
                Nperm mm3 iqr3(:)' 3;
                Nperm mm4 iqr4(:)' 4];
    dat = [dat;datloc];
end

writematrix(dat,sprintf('../data/compTimes_varyPerm_method%s.csv',tmeth));

function Jval = Jfun(tt,freq,acro,Qform,Amp,Nsamp)
	tt    = reshape(tt,[],1);

	mu    = Amp*cos(2*pi*freq.*tt -acro);
	x     = mu + randn([length(tt),1,1,1,Nsamp]);
	[pwr2,sgn] = evalChebPowerbnd(Qform,x,0.05,'rig');
	pwr2 = pwr2.*((-1).^(~sgn));
	Jval   =min(pwr2,[],'all');
end

function pwr = estTinfloc(tt,Nsamp,freq,acro,Amp,Nperm,Qform)
mu      = Amp*cos(2*pi*freq*tt -acro); % simulated signal
mu      = reshape(mu,[],1);
x       = mu + randn([length(tt),1,1,1,Nsamp]);

alpha   = .05;
pwrGrid = squeeze(fastMCTinfpower(Qform,x,Nperm,alpha));
pwr     = min(pwrGrid,[],2);
end


function pwr = estT2loc(tt,Nsamp,freq,acro,Amp,Nperm,Qform)
mu      = Amp*cos(2*pi*freq*tt -acro); % simulated signal
mu      = reshape(mu,[],1);
x       = mu + randn([length(tt),1,1,1,Nsamp]);

alpha   = .05;
pwrGrid = squeeze(fastMCT2power(Qform,x,Nperm,alpha));
pwr     = min(pwrGrid,[],2);
end

function [avgTime, iqrTime] = repTimeFunc(funcHandle, numRuns, timeMethod, varargin)
    iqrTime = NaN(1,2);
    switch timeMethod
		case 'timeit'
			avgTime = timeit(@() funcHandle(varargin{:}));
        case 'tictoc'
            times = zeros(1, numRuns);
            fprintf('    on trial'),        
            for i = 1:numRuns
                fprintf(' %d',i)
                tStart = tic;
			    funcHandle(varargin{:});
			    times(i) = toc(tStart);
            end
            avgTime = mean(times);
            [~,iqrTime] = iqr(times);
    end
    fprintf('\n')
end
