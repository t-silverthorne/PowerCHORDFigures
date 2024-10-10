tic
Bmat_master=[];

windows = 12:2:24;
windows = windows-1;
nt      = 48;
N       = 8;
addpath('../../utils/')
    
tic
for window=windows
    window
    system('rm output_*');
    s1=strjoin({'awk -f "night_filt.awk"',num2str(window),'1'});
    s2=strjoin({'cNecks_48_',num2str(N),'.txt > sols_temp.txt'},'');
    system(strjoin({s1,s2}));    
    
    best_eig = -Inf;
    fileList = dir('output_*.txt');

    for ff=1:length(fileList)
        
        % read into matlab
        fileID       = fopen(fileList(ff).name, 'r');
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
        
        % find optimal
        [~,eig]  = getMinEigMulti(Smat',1,1,1,false);
        if max(eig)>best_eig
            [best_eig,mind] = max(eig);
            topt            = Smat(mind,:);
            bopt            = Bmat(mind,:);
        end
    end
    
    % compare to alternate practice of equispaced
    %polarplot(2*pi*topt,1,'.k')
    t_unif = ((1:N)/N-1/N)';
    fprintf("\nAlternative:  %f\n",getMinEig(t_unif,1));
    fprintf('Optimal    :  %f\n',max(eig));
    
    clear Bmat Smat Tmat 
    Bmat_master = [Bmat_master; bopt];
end
toc
Bmat_master

writematrix(Bmat_master,'../../../clean_figs/data/window_sols_fp2.csv')
toc
%example awk command: awk '{line=$0 $0; if (line ~ /000000000001000000000001/) print $0}' test_in.txt > test_out.txt