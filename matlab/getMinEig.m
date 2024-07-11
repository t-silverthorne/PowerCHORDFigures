function emax = getMinEig(t,freq,method)
% return min eigenvalue of B inverse
arguments
    t;
    freq;
    method='div';
end
A    = [0 0;1 0;0 1];
X    = [ones(length(t),1) cos(2*pi*freq*t) sin(2*pi*freq*t)];
B    = A'*((X'*X)\A);
switch method
    case 'div'
        emax = 1/max(eig(B));
    case 'inv'
        emax = min(eig(inv(B)));
    otherwise
        error('unknown method');
end
end

