function [Cm11,Cm12,Cm21,Cm22] = getFourQuadBlocks(n,Nm,freq)
% make the four block matrices that define quadratic program for YALMIP
tau  = ((1:n)/n -1/n)';
c    = reshape(cos(2*pi*freq*tau),[],1);
s    = reshape(sin(2*pi*freq*tau),[],1);

dCS  = [diag(c.*c) diag(c.*s);
        diag(c.*s) diag(s.*s) ];
CS   = [c*c' s*c'
        c*s' s*s'];
Cmat = (dCS-CS/Nm);

Cm11 = Cmat(1:n,1:n);
Cm12 = Cmat(1:n,(n+1):2*n);
Cm21 = Cmat((n+1):2*n,1:n);
Cm22 = Cmat((n+1):2*n,(n+1):2*n);
end

