clear all
n    = 144;
Nm   = 32;
tau  = ((1:n)/n -1/n)';
fvec = linspace(1,24,24);
sdpvar eta;
mu = binvar(n,1);

F=[sum(sum(mu))==Nm];

for freq=fvec
    c = reshape(cos(2*pi*freq*tau),[],1);
    s = reshape(sin(2*pi*freq*tau),[],1);
    
    dCS = [diag(c.*c) diag(c.*s);
           diag(c.*s) diag(s.*s) ];
    
    CS=[c*c' s*c'
        c*s' s*s'];
    
    IOmat = [ones(1,n) zeros(1,n);zeros(1,n) ones(1,n)];
    
    Cmat = (dCS-CS/Nm);
    Cm11 = Cmat(1:n,1:n);
    Cm12 = Cmat(1:n,(n+1):2*n);
    Cm21 = Cmat((n+1):2*n,1:n);
    Cm22 = Cmat((n+1):2*n,(n+1):2*n);

    F = [F, eta <=lambda_min([mu'*Cm11*mu mu'*Cm12*mu; mu'*Cm12*mu mu'*Cm22*mu])];
end
h =-eta;

options=sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpsolver','gurobi', ... % can be gurobi or mosek
                    'bmibnb.maxtime',10*60);
warning('off','all')
optimize(F,h,options)
warning('on','all')