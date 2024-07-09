addpath('../')
Nm     = 50;
fmin   = 1;
fmax   = 24;
Nbatch = 100;
t    = rand(Nm,1,1,Nbatch);
tic;getMinEigMulti(t,fmin,fmax,1024,true);toc
tic;getMinEigMulti(t,fmin,fmax,1024,false);toc


