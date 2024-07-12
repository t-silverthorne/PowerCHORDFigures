Nm    = 48;
Npop  = 1e2;
fmin  = 1;
fmax  = 24;
Nfreq = 2^10;

eps = .05;
Tmat     = rand(Nm,Npop);

useGPUglob=true;
Niter = 500;
tic;
scores=NaN(1,Niter)


for ii=1:Niter
    % score population
    [~,Jnow] = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
    scores(ii) = max(Jnow);
    % select parents
    cind      = cell2mat(arrayfun(@(ii) randsample([1:(ii-1) (ii+1):Npop], ...
                                       3,true,Jnow([1:(ii-1) (ii+1):Npop]))', ...
                         1:Npop,'UniformOutput',false));
    % score candidates
    Tcand     = Tmat(:,cind(1,:)) + eps*(Tmat(:,cind(2,:)) - Tmat(:,cind(3,:)));
    Tcand     = mod(Tcand,1);
    [~,Jcand] = getMinEigMulti(Tcand,fmin,fmax,Nfreq,useGPUglob);
    % evolve population
    Tmat(:,Jcand>Jnow) = Tcand(:,Jcand>Jnow);
end
toc
plot(1:Niter,scores,'.k')