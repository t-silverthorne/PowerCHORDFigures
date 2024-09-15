directory = 'sweep_yalmip_etabd';
matFiles = dir(fullfile(directory, '*.mat'));
n=144;
nodes=[]
for ii = 1:length(matFiles)
    % Load the .mat file
    filePath = fullfile(directory, matFiles(ii).name);
    data = load(filePath);
    nodes(end+1) = data.res{4}.solveroutput.nodes;
    if data.res{4}.solveroutput.nodes==1
        data.res{4}.solveroutput.lower == -data.res{6} 
    end
end