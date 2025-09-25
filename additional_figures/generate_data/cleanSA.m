%% ---------- parameters ---------------
% clean wrapper for simulated annealing
tic;
clear;rng('default');
clf;
addpath('../../MATLAB/utils')
Nmeas   = 48; % number of measurements
Amp     = 1;
MaxIter = 10000;
fmin    = 1;         % min freq in window
fmax    = Nmeas/2;   % max freq in window
tt      = linspace(0,1,Nmeas+1);
tt      = tt(1:end-1)';
mode    = 'test';
switch mode % 24 5 test shows repetition in measurement pattern
    case 'test'
        % cheb params
        Nfreq_ch = 48;   % num freqs for Cheb bound
        Nacro_ch = 8;   % num acros for Cheb bound
        Nsamp_ch = 5e1; % for Cheb bound
        Nfq_Tinf = 8;   % num freqs for constructing test statistic
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nsamp_mc = 1e2; 
        Nfreq_mc = 48;
        Nacro_mc = 8;
        Nperm_mc = 1e2; 
    case 'test2'
        % cheb params
        Nfreq_ch = 8;   % num freqs for Cheb bound
        Nacro_ch = 8;   % num acros for Cheb bound
        Nsamp_ch = 1e2; % for Cheb bound
        Nfq_Tinf = 8;   % num freqs for constructing test statistic
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nsamp_mc = 1e3; 
        Nfreq_mc = 16;
        Nacro_mc = 8;
        Nperm_mc = 5e2; 


    case 'moderate'
        % cheb params
        Nfreq_ch = 32;  % num freqs for Cheb bound
        Nacro_ch = 16;  % num acros for Cheb bound
        Nsamp_ch = 5e1; % for Cheb bound
        Nfq_Tinf = 32;  % num freqs for constructing test statistic
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nsamp_mc = 1e2; 
        Nfreq_mc = 16;
        Nacro_mc = 16;
        Nperm_mc = 1e2;

	case 'real'
        % cheb params
        Nfreq_ch = 64;  % num freqs for Cheb bound
        Nacro_ch = 64;  % num acros for Cheb bound
        Nsamp_ch = 5e1; % for Cheb bound
        
		Nfq_Tinf = 64;  % num freqs for constructing test statistic
        Nfq_T2   = 1e3; % num freqs for constructing test statistic
        
        % Monte Carlo params
        Nsamp_mc = 2e2; 
        Nfreq_mc = 64;
        Nacro_mc = 64;
        Nperm_mc = 2e2; 
end
% ---------- pre optimization ----------------
tiledlayout(2,2);

[pwr2_mc,pwrinf_mc,pwr2_ch,fmc,fch]=benchmarkDesign(tt,fmin,fmax,Amp,...
                Nfreq_ch,Nacro_ch,Nsamp_ch,Nfq_Tinf,Nfq_T2, ...
                Nfreq_mc,Nacro_mc,Nsamp_mc,Nperm_mc,fmin,fmax*.95);
fprintf('Equispaced power Tinf     MC:   %d\n',min(pwrinf_mc))
fprintf('Equispaced power T2       MC:   %d\n',min(pwr2_mc))
fprintf('Equispaced power T2       CH:   %d\n',min(pwr2_ch))
nexttile(1)
plot(tt,1,'.k')
nexttile(3)
fprintf("----\n")
% plot equispaced
plot(fmc,pwrinf_mc,'-k')
hold on
plot(fmc,pwr2_mc,'-b')
plot(fch,pwr2_ch,'--b')

% ---------- optimization --------------------
freqs_ch = linspace(fmin,fmax,Nfreq_ch);
acros_ch = rand(Nacro_ch,1)*2*pi;
freqs_ch = reshape(freqs_ch,1,1,1,1,1,[]);
acros_ch = reshape(acros_ch,1,1,1,1,1,1,[]);
fqf_2    = reshape(linspace(fmin,fmax,Nfq_T2),1,1,[]);

Jwrap = @(tt) - Jfun(tt,freqs_ch,acros_ch,fqf_2,Amp,Nsamp_ch);
Jwrap(tt)

tt0 = linspace(0,1,Nmeas+1);
tt0 = tt0(1:Nmeas)';
lb = zeros(size(tt0));
ub = ones(size(tt0));
options = optimoptions('simulannealbnd', ...
    'Display','iter', ...
    'MaxIterations',MaxIter, ...
    'MaxFunctionEvaluations',1e4);

[tt_opt,fval,exitflag,output] = simulannealbnd(Jwrap,tt0,lb,ub,options)

% ---------- post optimization ---------------

fprintf('Equispaced power Tinf     MC:   %d\n',min(pwrinf_mc))
fprintf('Equispaced power T2       MC:   %d\n',min(pwr2_mc))
fprintf('Equispaced power T2       CH:   %d\n',min(pwr2_ch))
nexttile(2)
tt = tt_opt;
plot(tt,1,'.k')
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
	%xlim([0,1]);
	ylim([0,1.1]);
end
toc
%fname = sprintf('test_Nmeas%d_Amp%d_MaxIter%d_fmin%d_fmax%d_Nfreqch%d_Nacroch%d_Nsampch%d_NfqTinf%d_NfqT2%d_Nsampmc%d_Nfreqmc%d_Nacromc%d_Npermmc%d.fig', ...
%    Nmeas, Amp, MaxIter, fmin, fmax, ...
%    Nfreq_ch, Nacro_ch, Nsamp_ch, Nfq_Tinf, Nfq_T2, ...
%    Nsamp_mc, Nfreq_mc, Nacro_mc, Nperm_mc);
%savefig(fname);
function Jval = Jfun(tt,freqs,acros,fqf_2,Amp,Nsamp)
	tt    = reshape(tt,[],1);
	Q2    = getQuadForm(tt,fqf_2);

	mu    = Amp*cos(2*pi*freqs.*tt -acros);
	sz    = size(mu);
	x     = mu + randn([sz(1:4),Nsamp,sz(6:end)]);
	[pwr2,sgn] = evalChebPowerbnd(Q2,x,0.05,'rig');
	pwr2 = pwr2.*((-1).^(~sgn));
	Jval   =min(pwr2,[],'all');
end
