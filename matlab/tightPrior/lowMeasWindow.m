addpath('../')
fmin  = 1;
fmax  = 1;
Nfreq = 1;

Nm    = 6; % number during the day

settings.M          = 0; % number of measurements during the night
settings.Npop       = 1e2;
settings.Niter      = 1e3;
settings.eps        = .05;
settings.useGPUglob = true;
settings.time_max   = Inf;
settings.CR         = .05;

M=settings.M;
t1      = (1:Nm)/Nm - 1/Nm;
t2      = (1:M)/M-1/M;
t1      = reshape(t1,[],1);
t2      = 0.5+0.5*reshape(t2,[],1);
[~,Jnow]         = getMinEigMulti(toState(t1,t2),fmin,fmax,Nfreq,true)

[Tmat,eigfinal,scores]=pCHORDlowmeas(Nm,fmin,fmax,Nfreq,settings);

tiledlayout(2,1)
nexttile
plot(scores,'.k')
yline(Jnow,'-k')
ylim([0,10])
nexttile
[~,ind]=max(eigfinal);
plot(toState(Tmat(:,ind),t2),1,'.k')
hold on
plot(toState(((1:Nm)/Nm-1/Nm)',t2),2,'.k')
xlim([0,1])
ylim([0,3])
%%
toState(Tmat(:,ind),t2)
toState(((1:Nm)/Nm-1/Nm)',t2)
%%
function Smat = toState(Tmat,t2)
Smat=vertcat(0.5*Tmat,repmat(t2,1,size(Tmat,2)));
end

function [Tmat,eigfinal,scores] = pCHORDlowmeas(Nm,fmin,fmax,Nfreq,settings)

    % unpack
    Npop       = settings.Npop;
    Niter      = settings.Niter;
    time_max   = settings.time_max;
    eps        = settings.eps;
    CR         = settings.CR;
    M          = settings.M;
    useGPUglob = settings.useGPUglob;

    t2      = (1:M)/M-1/M;
    t2      = 0.5*t2 + 0.5;
    t2      = reshape(t2,[],1);
    
    Tmat       = rand(Nm,Npop);
    scores     = [];
    tic
    ii=1;
    while (ii <=Niter) && (toc<time_max)
        ii
        % score population
        [~,Jnow]         = getMinEigMulti(toState(Tmat,t2),fmin,fmax,Nfreq,useGPUglob);
        scores(ii)       = max(Jnow);
        
        % evolution
        cind             = cell2mat(arrayfun(@(ii) randsample(1:Npop,3,false)',1:Npop,'UniformOutput',false));
        Tcand            = Tmat(:,cind(1,:)) + eps*(Tmat(:,cind(2,:)) - Tmat(:,cind(3,:)));
        Tcand            = mod(Tcand,1);

        % crossover
        Tcr              = Tmat;
        ind_cross        = rand(Nm,Npop)<CR;
        Tcr(ind_cross)   = Tcand(ind_cross);
        
        
        % evolve population
        [~,Jcr]          = getMinEigMulti(toState(Tcr,t2),fmin,fmax,Nfreq,useGPUglob);
        Tmat(:,Jcr>Jnow) = Tcr(:,Jcr>Jnow);

        ii=ii+1;
    end
    [~,eigfinal] = getMinEigMulti(toState(Tmat,t2),fmin,fmax,Nfreq,useGPUglob);
end