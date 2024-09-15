N = 6;
nrep =1e5;
uvec=NaN(1,nrep);
for ii=1:nrep
    
    
    b = rand(2,1);
    D = rand(2,2);
    D = D*D';
    M1 = [N b'
         b D];
    M2 = [1 b'/sqrt(N)
         b/sqrt(N) D];
    logic1=min(eig(M1))>0;
    logic2=min(eig(M2))>0;

    uvec(ii)=logic1==logic2;
end

prod(uvec)


