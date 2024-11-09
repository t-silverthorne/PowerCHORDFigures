%% Extract results of differential evolution

directory = '../../PowerCHORD/MATLAB/cluster_compute_gpu/diffEvolveOutput/';
matFiles = dir(fullfile(directory, '*.mat'));
%%
T=table();
prob_inds=[]
for ii = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(directory, matFiles(ii).name);
    try
        data = load(filePath);
    catch
        prob_inds(end+1)=ii;
    end
    
    sp      = strsplit(matFiles(ii).name,'_');
    Nmeas   = str2num(sp{2});
    fmin    = str2num(sp{4});
    fmax    = str2num(sp{6});
    method  = "diffEVCR";
    [ncp,ind] = max(data.eigfinal);
    tvec    = data.Tmat(:,ind)';
    upper   = NaN;
    MIPgap  = NaN;
    lpred   = NaN;
    tvec    = [tvec NaN(1,48-length(tvec))];
    T       =  [T;table(method,Nmeas,fmin,fmax,ncp,upper,lpred,MIPgap,tvec)];
end
prob_inds
%%
writetable(T,'../data/diffEvolveOutput.csv')
