function filterNecklaces()
%
global Nmat Bmat

n     = size(Nmat,2);
rNmat = Nmat(:,n:-1:1);

Bmat=Nmat;
ii=1;
while ii <=size(rNmat,1)
    rNmat(ii,:)=minRotation(rNmat(ii,:));
    if any(rNmat(ii,:)~=Bmat(ii,:))
        idx=all(Bmat==rNmat(ii,:),2);
        Bmat(idx,:)=[];
        rNmat(idx,:)=[];
    end
    ii=ii+1;
end
end

