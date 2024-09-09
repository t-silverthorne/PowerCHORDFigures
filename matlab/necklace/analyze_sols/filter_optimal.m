function bopt = filter_optimal(N,m,nt,solsdir)
% filter for measurement schedules containing a specific substring
%  calculate power for that family of measurement schedules
%  return optimal 

    % construct filename
    file_in  = strcat('cNecks_',string(nt),'_',string(N+m),'.txt');
    file_in  = strcat(solsdir,file_in);
    file_out = strjoin(['temp',string(N),string(m),string(nt),'.txt'],'_');
    
    % construct query string
    qt       = zeros(1,nt/2); 
    qi       = linspace(1,nt/2+1,m+1);
    qi       = qi(1:end-1);
    qt(qi)   = 1;
    qt       = strjoin(string(qt),'');
    qt       = strcat(' /',qt,'/');
    
    % query solution database using awk command 
    awk_1    = "awk '{line=$0 $0; if (line ~ ";
    awk_2    = strcat(qt,") print $0}'");
    awk_3    = strcat(file_in," > ", file_out);
    awk_cmd  = strcat(awk_1,awk_2,awk_3);
    awk_cmd2 = strcat(awk_1,awk_2,file_in);
    system(awk_cmd);
    system(strjoin(["wc -l ",file_out],''));
    
    % read into matlab
    fileID       = fopen(file_out, 'r');
    data         = textscan(fileID, '%s');
    fclose(fileID);
    charData     = char(data{1});
    Bmat         = charData - '0';  % Subtract '0' to convert char '0'/'1' to numeric 0/1
    
    % construct time matrix
    tvec = (1:nt)/nt - 1/nt;
    Tmat = repmat(tvec,size(Bmat,1),1);
    Tmat = Tmat.*Bmat;
    Smat = zeros(size(Tmat,1),N+m);
    for ii=1:size(Smat,1)
        Smat(ii,:) = Tmat(ii,find(Tmat(ii,:)));
    end
    
    % find optimal
    [~,eig]  = getMinEigMulti(Smat',1,1,1,false);
    [~,mind] = max(eig);
    bopt     = Bmat(mind,:);



end

