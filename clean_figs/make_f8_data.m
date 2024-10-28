nvec = [24,48,96]
cMat = NaN(3,17);

for ii=1:length(nvec)
    cMat(ii,:)=arrayfun(@(dd) enum(nvec(ii),dd),4:20)
end

writematrix(log10(cMat),'data/neck_counts.csv')

%%
plot(4:20,log10(cMat(1,:)))



