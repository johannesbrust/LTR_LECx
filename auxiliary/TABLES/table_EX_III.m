%-------------------- table_EX_III ---------------------------------------%
%
% Script to store data from Experiment into a table for the
% manuscript
%
% The table has 50 rows and 2 + 2 * 7 = 16 columns
%
%-------------------------------------------------------------------------%
% 10/16/20, J.B.
% 10/21/20, J.B., problem mask
% 06/07/21, J.B., Modification highlighting results
% 06/09/21, J.B., Modification fo release

% Load data from experiment
%datapath = fullfile(pwd,'..','..','/data/EXPERIMENT_III_EXT');
datapath = fullfile(pwd,'..','..','/data/');
dataEX = load([datapath,'EXPERIMENT_III_EXT']);

% Load problem data
dataProb = load('EVEN_N_PROBS');
probMask = dataProb.probMask;

nrow    = 50; % 30
%ncol    = 4 + 3*9;
table = cell(nrow,1);
nsol    = 7;
skipIdx = 0; % Don't include reprojected solver

for i = 1:nrow
    
    % Process data
    pname_L = dataProb.names{probMask(i)};
    pname_SP= strsplit(pname_L,'/');
    pname = pname_SP{2};
    
    pname = strrep(pname,'_','\_');
    m = dataProb.ms(probMask(i),1);
    n = dataProb.ms(probMask(i),2);
    rnk = dataProb.msRks(probMask(i),2);
    dens = (dataProb.nnzs(probMask(i))/(m*n));
    
    tblr = sprintf('$\\texttt{%s}$ & %i/%i & %i/%.1g &',pname,m,n,rnk,dens);
    %dataEX.conds(i));
    
    % Sort outcomes
    %sit = sort(dataEX.its(i,:),'ascend');
    sti = sort(dataEX.times(i,:),'ascend');
    
    for j = 1:nsol
        
        lastc = '&';
        if j == nsol
            lastc ='\\';
        end
        
        if j~=skipIdx
            % Different row depending on outcome
            ex = dataEX.exs(i,j); % Either 0 or 1
            conv = dataEX.convs(i,j);
            if ex == 0 && isnan(conv) == 1
                tblr = [tblr, sprintf('$\\texttt{N/A}^{*}$ & $\\texttt{N/A}$ %s',lastc)]; %#ok<AGROW>
            elseif ex == 0 && isnan(conv) == 0
                tblr = [tblr, sprintf('$\\texttt{NC}^{\\dagger}$ & $\\texttt{NC}$ %s',lastc)]; %#ok<AGROW>
            else
                it = dataEX.its(i,j);
                ti = dataEX.times(i,j);
                tfrmt = '%.2g';
                %                         if ti > 1e2;
                %                             tfrmt = '%.0f';
                %                         end
                
                if ti == sti(1)
                    tfrmt = '\\textbf{%.2g}';
                    if ti > 1e2;
                        tfrmt = '\\textbf{%.0f}';
                    end
                elseif ti == sti(2)
                    tfrmt = '\\emph{%.2g}';
                    if ti > 1e2;
                        tfrmt = '\\emph{%.0f}';
                    end
                end
                
                tblr = [tblr, sprintf(['%i &',tfrmt,'%s'],...
                    round(it),ti,lastc)];  %#ok<AGROW>
                %                 tblr = [tblr, sprintf('%i & %.2g %s',...
                %                         round(dataEX.its(i,j)),dataEX.times(i,j),lastc)];  %#ok<AGROW>
            end
            
        end
        
    end
    
    table{i} = tblr;%[tblr,'\\'];
    
end

% Write output table to file
fid_w   = fopen([datapath,'table_EX3.txt'],'w');
fprintf(fid_w,'%s \n',table{:});
fclose(fid_w);

% display(table)