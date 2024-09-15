function [prob,mu,eta,F] = pCHORD(n,Nm,fvec,options,method,settings)
arguments % can leave options blank and yet YALMIP decide on everything
    n;Nm;fvec;
    options  = sdpsettings();
    method   = 'YALMIP';
    settings = struct();
end



% defaults 
if isempty(fieldnames(settings))
    switch method
        case 'YALMIP'
            settings.warmstart         = 'auto'; % options: none, user, auto
            settings.relax             = false;
            settings.fix_gauge         = true;
            settings.force_meas_region = false;
            settings.no_meas_region    = false;
            settings.FMRvec            = zeros(n,1);
            settings.NMRvec            = zeros(n,1);
        otherwise
            error('Method not implemented')
    end
end

switch method
    case 'YALMIP'
        if settings.fix_gauge && (settings.force_meas_region || settings.no_meas_region)
            error('Cannot fix gauge and enforce mandatory/no meas region')
        end
        % decision variables
        if settings.relax
            mu = sdpvar(n,1);
            F = [mu>=0,mu<=1];
        else
            mu = binvar(n,1);
            F  = [];
        end
        eta=sdpvar(1);
        
        % constraints
        F  = [F,sum(sum(mu))==Nm];
   	    F  = [F,eta<=Nm/2,0<=eta];

        for freq=fvec
            [Cm11,Cm12,Cm21,Cm22] = getFourQuadBlocks(n,Nm,freq);
            F = [F, eta <=lambda_min([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm21*mu mu'*Cm22*mu])];
        end
        if settings.fix_gauge
            th1=(1:n)/n - 1/n;
            F=[F,-1/n <= sum(diag(sin(2*pi*th1))*mu)/n <= 1/n];
        end
        if settings.force_meas_region
            F=[F,settings.FMRvec'*mu == sum(settings.FMRvec)];
        end

        if settings.no_meas_region
            F=[F,settings.NMRvec'*mu == 0];
        end

        % objective
        h  = -eta; 

        switch settings.warmstart
            case 'auto'
                mu0 = zeros(n,1);
                mu0(randsample(1:n,Nm,false))=1;  
            case 'user'
                mu0=settings.mu0;
        end
        switch settings.warmstart
            case {'auto','user'}
               if settings.fix_gauge
                    th1=(1:n)/n-1/n;
                    for ii=1:n
                        mu0 = [mu0(end); mu0(1:end-1)];
                        sh_list{ii}=mu0;
                        score(ii) = sin(2*pi*th1)*mu0/n;
                    end
                    [~,ind]=min(abs(score));
                    mu0 = sh_list{ind};
                end
                eta0=[];
                for freq=fvec
                    [Cm11,Cm12,Cm21,Cm22] = getFourQuadBlocks(n,Nm,freq);
                    eta0(end+1) = min(eig([mu0'*Cm11*mu0 mu0'*Cm12*mu0; mu0'*Cm21*mu0 mu0'*Cm22*mu0]));
                end
                eta0=min(eta0);
                
                warmstart(eta,eta0);
                warmstart(mu,mu0);
        end

        warning('off','all') % TODO: make this optional in future
        prob=optimize(F,h,options);% solve
        warning('on','all')

    otherwise
        error('Method not implemented')
end


end

