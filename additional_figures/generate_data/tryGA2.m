%% ---------- parameters ---------------
% use genetic algorithm to enforce linear constraints
% delete(gcp('nocreate')) % don't need to run every time
% parpool('local',6);  

tic;
clear;rng('default');
clf;
addpath('../MATLAB/utils')
PopSize = 100;
dx      = .1;
Nmeas   = 48; % number of measurements
ceps    = .05/Nmeas;
Amp     = 2;
MaxIter = 10;
fmin    = 1;         % min freq in window
fmax    = Nmeas/2;   % max freq in window
mode    = 'tiny';
switch mode % 24 5 test shows repetition in measurement pattern
    case 'tiny'
        % cheb params
        Nfreq_ch = 48;   % num freqs for Cheb bound
        Nacro_ch = 8;   % num acros for Cheb bound
        Nsamp_ch = 1e1; % for Cheb bound
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nfq_Tinf = 48;   % num freqs for constructing test statistic
        Nsamp_mc = 1e2; 
        Nfreq_mc = 48;
        Nacro_mc = 8;
        Nperm_mc = 20; 
    case 'test'
        % cheb params
        Nfreq_ch = 48;   % num freqs for Cheb bound
        Nacro_ch = 16;   % num acros for Cheb bound
        Nsamp_ch = 1e1; % for Cheb bound
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nfq_Tinf = Nmeas;   % num freqs for constructing test statistic
        Nsamp_mc = 1e2; 
        Nfreq_mc = 48;
        Nacro_mc = 8;
        Nperm_mc = 50; 
	case 'real'
        % cheb params
        Nfreq_ch = 64;  % num freqs for Cheb bound
        Nacro_ch = 64;  % num acros for Cheb bound
        Nsamp_ch = 1e1; % for Cheb bound
        
		Nfq_Tinf = 128;  % num freqs for constructing test statistic
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nsamp_mc = 2e2; 
        Nfreq_mc = 64;       
        Nacro_mc = 64;
        Nperm_mc = 2e2; 
end

freqs_ch = linspace(fmin,fmax,Nfreq_ch);
acros_ch = rand(Nacro_ch,1)*2*pi;
freqs_ch = reshape(freqs_ch,1,1,1,1,1,[]);
acros_ch = reshape(acros_ch,1,1,1,1,1,1,[]);
fqf_2    = reshape(linspace(fmin,fmax,Nfq_T2),1,1,[]);


initPop = NaN(PopSize,Nmeas);
for jj =1:PopSize   
    initPop(jj,:) = reshape(randInitDesign(Nmeas,ceps),1,[]);
end

Acstr = orderConstraintMat(Nmeas);
bsctr = -ceps*ones(Nmeas-1,1);

Aeq = [1 zeros(1,Nmeas-1); zeros(1,Nmeas/2) 1 zeros(1,Nmeas/2-1)];
beq = [0;.5];
% Aeq = zeros(4,Nmeas);
% idxs = Nmeas/4 * (0:3) + 1;   % positions: 1, Nmeas/4+1, Nmeas/2+1, 3Nmeas/4+1
% for k = 1:4
%     Aeq(k,idxs(k)) = 1;
% end
% Aeq
% beq = [0;.25;.5;.75];

lb = zeros(Nmeas,1);
ub = ones(Nmeas,1);
Jwrap = @(tt) - Jfun(tt,freqs_ch,acros_ch,fqf_2,Amp,Nsamp_ch);

% GA options
options = optimoptions('ga', ...
    'Display','iter', ...
    'MaxGenerations',MaxIter, ...
    'MaxStallGenerations',50, ...
    'PopulationSize',PopSize, ...     
    'FunctionTolerance',1e-6, ...
    'UseParallel',false);
% 'InitialPopulationMatrix', initPop,...

% Run GA: need number of variables = length(tt0)
nvars = Nmeas;
tt_opt = ga(Jwrap,Nmeas,Acstr,bsctr,Aeq,beq,lb,ub,[],options)

% ---------- post optimization ---------------

% fprintf('Equispaced power Tinf     MC:   %d\n',min(pwrinf_mc))
% fprintf('Equispaced power T2       MC:   %d\n',min(pwr2_mc))
% fprintf('Equispaced power T2       CH:   %d\n',min(pwr2_ch))

nexttile(2)
tt = tt_opt';
plot(tt,1,'.k')
xlim([0,1])
nexttile(4)
[pwr2_mc,pwrinf_mc,pwr2_ch,fmc,fch]=benchmarkDesign(tt,fmin,fmax,Amp,...
                Nfreq_ch,Nacro_ch,Nsamp_ch,Nfq_Tinf,Nfq_T2, ...
                Nfreq_mc,Nacro_mc,Nsamp_mc,Nperm_mc,fmin,fmax);
fprintf('Optimal design power Tinf MC:   %d\n',min(pwrinf_mc))
fprintf('Optimal design power T2   MC:   %d\n',min(pwr2_mc))
fprintf('Optimal design power T2   CH:   %d\n',min(pwr2_ch))

plot(fmc,pwrinf_mc,'-k')
hold on
plot(fmc,pwr2_mc,'-b')
plot(fch,pwr2_ch,'--b')

for jj=1:4
	nexttile(jj)

	ylim([0,1.1]);
end
toc

function Jval = Jfun(tt,freqs,acros,fqf_2,Amp,Nsamp)
	tt    = reshape(tt,[],1);
	Q2    = getQuadForm(tt,fqf_2);

	mu    = Amp*cos(2*pi*freqs.*tt -acros);
	sz    = size(mu);
	x     = mu + randn([sz(1:4),Nsamp,sz(6:end)]);
	[pwr2,sgn] = evalChebPowerbnd(Q2,x,0.05,'rig');
	pwr2 = pwr2.*(sgn);
    Jval = min(pwr2,[],'all')+1e-3*mean(pwr2,'all');
end
