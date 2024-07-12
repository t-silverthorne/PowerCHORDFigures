function [prob,mu,eta] = pCHORD(n,Nm,fvec,options,relax)
arguments % can leave options blank and yet YALMIP decide on everything
    n;Nm;fvec;
    options = sdpsettings();
    relax   = false;
end

sdpvar eta; % decision variables
if relax
    mu = sdpvar(n,1);
else
    mu = binvar(n,1);
end
   
mu0 = zeros(n,1);
mu0(randsample(1:n,Nm,false))=1;

F  = sum(sum(mu))==Nm; % constraints

if relax
    F = [F,mu>=0,mu<=1]
end

eta0=[];
for freq=fvec
    [Cm11,Cm12,Cm21,Cm22] = getFourQuadBlocks(n,Nm,freq);
    F = [F, eta <=lambda_min([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm21*mu mu'*Cm22*mu])];
    eta0(end+1) = min(eig([mu0'*Cm11*mu0 mu0'*Cm12*mu0; mu0'*Cm21*mu0 mu0'*Cm22*mu0]));
end
eta0=min(eta0);

h =-eta; % objective to be minimized (YALMIP only does minimzation)

warmstart(mu,mu0);
warmstart(eta,eta0);

warning('off','all') % TODO: make this optional in future
prob=optimize(F,h,options);% solve
warning('on','all')
end

