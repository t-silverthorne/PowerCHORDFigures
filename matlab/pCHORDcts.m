function [Tmat,eigfinal,scores] = pCHORDcts(Nm,fmin,fmax,Nfreq,settings)

switch settings.method
    case 'diffEvolve'
        % unpack
        Npop       = settings.Npop;
        Niter      = settings.Niter;
        eps        = settings.eps;
        useGPUglob = settings.useGPUglob;
        
        % initial population
        Tmat       = rand(Nm,Npop);
        scores     = NaN(1,Niter);
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
        [~,eigfinal] = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
        
    otherwise
        error('method not yet implemented')
end
end

