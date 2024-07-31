clear

options=sdpsettings('solver','bmibnb', ...
                    'bmibnb.uppersolver','fmincon', ...   % can't be mosek
                    'bmibnb.lowersolver','mosek', ...   % has to be mosek
                    'bmibnb.lpreduce',0,...
                    'bmibnb.maxtime',60, ...
                    'usex0',1);

fvec = 1:.5:12;
nf = length(fvec);
n  = 144;
m  = 24;
mu = binvar(n,1);
sdpvar eta

th1=(1:n)/n - 1/n;
obj = -eta;

init_state='rand';
switch init_state
    case 'bunched'
        mu0 = [ones(m/2,1); zeros(n-m,1); ones(m/2,1)];
    case 'rand'
        mu0 = zeros(n,1);
        mu0(randsample(1:n,m,false))=1;
        for ii=1:n
            mu0 = [mu0(end); mu0(1:end-1)];
            sh_list{ii}=mu0;
            score(ii) = sin(2*pi*th1)*mu0/n;
        end
        [mi,ind]=min(abs(score));
        mu0 = sh_list{ind};
        %abs(sin(2*pi*th1)*mu0/n)
        %mi
end







F1 = sum(mu)==m;
M = [];
eta0 = [];


for freq=fvec
    [C11,C12,C21,C22] = getFourQuadBlocks(n,m,freq);
    A=[mu'*C11*mu mu'*C12*mu; mu'*C21*mu mu'*C22*mu];
    
    F1=[F1,A-eta*eye(2)>=0];
    M = blkdiag(M,A);
        
    eta0(end+1) = min(eig([mu0'*C11*mu0 mu0'*C12*mu0; mu0'*C21*mu0 mu0'*C22*mu0]));
end

F0 = F1;
F1 = [F1,-2/n <= sum(diag(sin(2*pi*th1))*mu)/n <= 2/n];
F2 = [sum(mu)==m,M-eta*eye(2*nf)>=0];
F2 = [F2,-2/n <= sum(diag(sin(2*pi*th1))*mu)/n <= 2/n];

assign(mu,mu0);
assign(eta,min(eta0));
prob0=optimize(F0,obj,options);
v0  =value(obj);

assign(mu,mu0);
assign(eta,min(eta0));% todo check that it complains if you supply nonsense
prob1=optimize(F1,obj,options);
v1 = value(obj);


assign(mu,mu0);
assign(eta,min(eta0));
prob2=optimize(F2,obj,options);
v2=value(obj);

v0
v1
v2
prob0
prob1
prob2
%%

