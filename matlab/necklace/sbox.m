'starting'
clear all
global Nmat
Nmat=[];
a(1)=0;

tic
gen(1,1,4,a);
toc

Nmat

%%
clear all
global Nmat
Nmat=[];
a(1)=1;
b(1)=1;
genfix(1,1,7,3,a,b)