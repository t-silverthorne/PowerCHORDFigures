%% Diffevolve sweep

directory = 'sweep_diffevolve';
matFiles = dir(fullfile(directory, '*.mat'));

T=table();
% Loop through each .mat file
for ii = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(directory, matFiles(ii).name);
    data = load(filePath);
    
    sp      = strsplit(matFiles(ii).name,'_');
    Nmeas   = str2num(sp{2});
    fmin    = str2num(sp{4});
    fmax    = str2num(sp{6});
    method  = "diffEV";
    [~,ind] = max(data.res{2});
    ncp     = data.res{2}(ind);
    tvec    = data.res{1}(:,ind)';
    upper   = NaN;
    MIPgap  = NaN;
    lpred   = NaN;
    tvec    = [tvec NaN(1,48-length(tvec))];
    T       =  [T;table(method,Nmeas,fmin,fmax,ncp,upper,lpred,MIPgap,tvec)];
end

% Yalmip

directory = 'yalmip_sweep';
matFiles = dir(fullfile(directory, '*.mat'));
n=144;
for ii = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(directory, matFiles(ii).name);
    data = load(filePath);
    
    sp      = strsplit(matFiles(ii).name,'_');
    Nmeas   = str2num(sp{2});
    fmin    = str2num(sp{4});
    fmax    = str2num(sp{6});
    method  = "YALMIP";
    ncp     = data.res{6};
    upper   =  - data.res{4}.solveroutput.lower;
    MIPgap  = 100*abs(upper-ncp)/abs(ncp);
    mu      = data.res{5}';
    tau     = (1:n)/n -1/n;
    tvec    = tau(mu>0);
    lpred   = ~sum(contains(filePath,'lpred_off'));
    tvec    = [tvec NaN(1,48-length(tvec))];
    T       =  [T;table(method,Nmeas,fmin,fmax,ncp,upper,lpred,MIPgap,tvec)];
end
writetable(T,'../solutions/powerCHORD_sols.csv')
