function [mu,eta] = pCHORD(n,Nm,fvec,options)
arguments % can leave options blank and yet YALMIP decide on everything
    n;Nm;fvec;
    options=sdpsettings();
end

sdpvar eta; % decision variables
mu = binvar(n,1);

F  = sum(sum(mu))==Nm; % constraints
for freq=fvec
    [Cm11,Cm12,Cm21,Cm22] = getFourQuadBlocks(n,Nm,freq);
    F = [F, eta <=lambda_min([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm21*mu mu'*Cm22*mu])];
end

h =-eta; % objective to be minimized (YALMIP only does minimzation)

warning('off','all') % TODO: make this optional in future
optimize(F,h,options) % solve
warning('on','all')
end

