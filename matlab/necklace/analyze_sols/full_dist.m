% construct filename
clear all
close all
solsdir  = ' ~/research/powerCHORD2/solutions/';
addpath('../../utils/')
system('rm *.txt')

nt   = 72;
freq = 1;
useGPU=false;
N    = 6;

file_in  = strcat('cNecks_',string(nt),'_',string(N),'.txt');
file_in  = strcat(solsdir,file_in);
file_out = strjoin(['tempDIST',string(N),string(nt)],'_');

fout=['"' 'output_' '"'];
dtxt=['"' '.txt' '"'];

len       =nt;
max_lines = 1e6;

% split the data
acmd = strjoin(["awk '", ...
"BEGIN {",  ...
"    file_num = 1", ...
"    line_count = 0", ...
strcat("    len =",string(len)) ,...
strcat("    max_lines =",string(max_lines)), ...
"}", ...
"length($0) == len {", ...
"    if (line_count >= max_lines) {", ...
"        file_num++", ...
"        line_count = 0", ...
"    }", ...
strcat("    print > " ,fout," file_num ",dtxt), ...
"    line_count++", ...
"}", ...
strcat("' ",file_in)],'\n');
system(acmd);

% load the data
eig_master = [];
split_files = dir('output_*.txt');
all_files   = {split_files.name};
tic
for file_loc=all_files
    fileID       = fopen(file_loc{1}, 'r');
    data         = textscan(fileID, '%s');
    fclose(fileID);
    charData     = char(data{1});
    Bmat         = charData - '0';  % Subtract '0' to convert char '0'/'1' to numeric 0/1
    
    % construct time matrix
    tvec = (1:nt)/nt - 1/nt;
    Tmat = repmat(tvec,size(Bmat,1),1);
    Tmat = Tmat.*Bmat;
    Smat = zeros(size(Tmat,1),N);
    for ii=1:size(Smat,1)
        Smat(ii,:) = Tmat(ii,find(Tmat(ii,:)));
    end
    
    % compute WCP
    [~,eig_loc]  = getMinEigMulti(Smat',freq,freq,1,useGPU);
    eig_master = [eig_master eig_loc];
end

size(eig_master)
tunif = (1:N)/N - 1/N;
[~,ref]=getMinEigMulti(tunif',freq,freq,1,false)
toc
histogram(eig_master)
xline(ref)

find(eig_master>ref+1e-12)