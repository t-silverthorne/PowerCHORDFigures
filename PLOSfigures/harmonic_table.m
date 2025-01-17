addpath('../PowerCHORD/MATLAB/utils/')
addpath('../PowerCHORD/MATLAB/optim_methods/')

%% find first harmonic which cannot be resolved
Nvec = 6:2:48;
nh   = 100;
hvec = NaN(1,length(Nvec));
for ii=1:length(Nvec)
    N = Nvec(ii);
    mt = (1:N)/N -1/N;
    mt = mt';
    [freqs,eigs]=getMinEigMulti(mt,1,nh,nh,false,'all','slow');
    bad_harm = find(abs((eigs-N/2)*2/N)>1e-13,1);
    if ~isempty(bad_harm)
        hvec(ii) = bad_harm;
    end
end
hvec
hvec

%% solve equivalent optimization problem and report power
ovec   = NaN(1,length(Nvec));
evec   = NaN(1,length(Nvec));
evec_s = NaN(1,length(Nvec));
evec_f = NaN(1,length(Nvec));
for ii=1:length(Nvec)
    if Nvec(ii)<12
        n=36;
    elseif Nvec(ii)<24
        n=48;
    else
        n=48*2;
    end
    [mu,eta,optim] = run_yalmip_timelim(1:hvec(ii),Nvec(ii),n,60*60);
    tau   = ((1:n)/n - 1/n)';
    ovec(ii)=optim.problem;
    evec(ii) = value(eta);
    [~,te]=getMinEigMulti(tau(value(mu)>.5),1,hvec(ii),hvec(ii),false,'min','slow');
    evec_s(ii) = te;
    [~,te]=getMinEigMulti(tau(value(mu)>.5),1,hvec(ii),hvec(ii),false,'min','fast');
    evec_f(ii) = te;
    save('table.mat')
		fid = fopen('progress.txt', 'w'); % Open file for writing
		fprintf(fid, 'Done iteration :\n');
		fprintf(fid, '%f\n', ii); % Print the variable elements one by one
		fclose(fid); % Close the file
end



%%

function [mu,eta,optim]=run_yalmip_timelim(freqs,Nmeas,n,tlim)
tau   = ((1:n)/n - 1/n)';

eta   = sdpvar(1,1,'full');
mu    = binvar(n,1);
F     = [sum(mu)==Nmeas,eta>=0,eta<=Nmeas/2];

for freq=freqs
    cvec = cos(2*pi*freq*tau);
    svec = sin(2*pi*freq*tau);
    uvec = ones(n,1);
    X=[uvec cvec svec];
    F=[F,X'*diag(mu)*X >= eta*eye(3)];
end
options=sdpsettings('solver','cutsdp');
options.cutsdp.maxtime=tlim;
optim=optimize(F,-eta,options);
check(F)
end
