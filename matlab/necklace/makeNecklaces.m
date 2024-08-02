function makeNecklaces(n,d)
%MAKENECKLACES Summary of this function goes here
%   Detailed explanation goes here

global Nmat
Nmat=[];
a=[0 0];

predictCount(n,d)

d=d+1;n=n+1;
jvals = ceil(n-d+1):-1:floor((n-1)/d+1);
tic
for j=jvals
    a(2) = j;
    b(2) = 1;
    genfix(1,1,n,d,a,b)
end
toc
size(Nmat)

fname = strcat('NeckMat_n_',string(n-1),'_d_',string(d-1));
save(fname,'Nmat');

end

