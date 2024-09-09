solsdir  = ' ~/research/powerCHORD2/solutions/';
nt=48;
smat = reshape(4:8,[],1);
smat = [smat repmat(2,length(smat),1) repmat(nt,length(smat),1)];
smat = [smat  NaN(size(smat,1),nt)];

for ii=1:size(smat,1)
    smat(ii,:) = [smat(ii,1:3) filter_optimal(smat(ii,1),smat(ii,2),smat(ii,3),solsdir)]
    pause(1)
end

writematrix(smat,strjoin(['optimal','nt',string(nt),'.csv'],'_'))

