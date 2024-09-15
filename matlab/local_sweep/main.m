
fmin = [1 2 6];
fmax = [12:4:24];
df   = .5;
n    = 144;
[fmin,fmax,Nmeas]=ndgrid(fmin,fmax,[16,32,48]);
tau   = ((1:n)/n - 1/n)';

pars=[fmin(:) fmax(:) Nmeas(:)];

size(pars)

ii      = 1

for ii=1:size(pars,1)
    fmin    = pars(ii,1);
    fmax    = pars(ii,2);
    Nmeas   = pars(ii,3);


    cvec = cos(2*pi*freq*tau);
    svec = sin(2*pi*freq*tau);

    

end