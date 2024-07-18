function [Tmat,eigfinal,scores] = pCHORDcts(Nm,fmin,fmax,Nfreq,settings)

switch settings.method
    case 'diffEvolveCR' % differential evolution with crossover
        
        % unpack
        Npop       = settings.Npop;
        Niter      = settings.Niter;
        time_max   = settings.time_max;
        eps        = settings.eps;
        CR         = settings.CR;
        useGPUglob = settings.useGPUglob;
        
        Tmat       = rand(Nm,Npop);
        scores     = [];
        tic
        ii=1;
        while (ii <=Niter) && (toc<time_max)
            % score population
            [~,Jnow]         = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
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
            [~,Jcr]          = getMinEigMulti(Tcr,fmin,fmax,Nfreq,useGPUglob);
            Tmat(:,Jcr>Jnow) = Tcr(:,Jcr>Jnow);

            ii=ii+1;
        end
        [~,eigfinal] = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
       

    case 'diffEvolve'
        % unpack
        Npop       = settings.Npop;
        Niter      = settings.Niter;
        time_max   = settings.time_max;
        eps        = settings.eps;
        useGPUglob = settings.useGPUglob;

        % initial population
        Tmat       = rand(Nm,Npop);
        scores     = [];
        tic
        ii=1;
        while (ii <=Niter) && (toc<time_max)
            % score population
            [~,Jnow] = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
            scores(ii) = max(Jnow);
            
            % select parents
            switch settings.diffEVsample
                case 'uniform'
                    cind      = cell2mat(arrayfun(@(ii) randsample(1:Npop, ...
                                               3,false)', ...
                                 1:Npop,'UniformOutput',false));                                    
                case 'fitness'
                    cind      = cell2mat(arrayfun(@(ii) randsample([1:(ii-1) (ii+1):Npop], ...
                                               3,true,Jnow([1:(ii-1) (ii+1):Npop]))', ...
                                 1:Npop,'UniformOutput',false));                    
            end

            % score candidates
            Tcand     = Tmat(:,cind(1,:)) + eps*(Tmat(:,cind(2,:)) - Tmat(:,cind(3,:)));
            Tcand     = mod(Tcand,1);
            [~,Jcr] = getMinEigMulti(Tcand,fmin,fmax,Nfreq,useGPUglob);
            
            % evolve population
            Tmat(:,Jcr>Jnow) = Tcand(:,Jcr>Jnow);

            ii=ii+1;
        end
        [~,eigfinal] = getMinEigMulti(Tmat,fmin,fmax,Nfreq,useGPUglob);
       
    case 'simulAnneal'
        % unpack
        useGPUglob = settings.useGPUglob;
        maxIter    = settings.Niter;
        maxTime    = settings.maxTime;
        cfun       = @(t) -Out2(@getMinEigMulti,t,fmin,fmax,Nfreq,useGPUglob,'min');
       
        % optimize
        opts       = optimoptions('simulannealbnd', ...
                     'MaxIterations',maxIter, ...
                     'HybridFcn','patternsearch','MaxTime',maxTime);
        t          = simulannealbnd(cfun,rand(Nm,1),zeros(Nm,1),ones(Nm,1),opts);
        
        % get output
        Tmat         = t;
        [~,eigfinal] = getMinEigMulti(t,fmin,fmax,Nfreq,useGPUglob,'min');
        scores       = []; %TODO: record evolution of eigenvalue during simulated annealing
    otherwise
        error('method not yet implemented')
end
end

