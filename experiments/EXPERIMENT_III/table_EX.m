%-------------------- table_EX -------------------------------------------%
%
% Script to store data from Experiment into a table for the
% manuscript
%
% The table has 30 rows and 3 + 2 * 7 = 17 columns
%
%-------------------------------------------------------------------------%
% 10/16/20, J.B.
% 10/21/20, J.B., problem mask

% Load data from experiment
datapath = fullfile(pwd,'..','..','/data/EXPERIMENT_III_EXT');
dataEX = load(datapath);

% Load problem data
dataProb = load('EVEN_N_PROBS');
probMask = dataProb.probMask;

nrow    = 50; % 30
%ncol    = 4 + 3*9;
table = cell(nrow,1);
nsol    = 7;
skipIdx = 0; % Don't include reprojected solver

% Data stored as (only rows 1-30 of the problems were used in experiment):
% dataEX = 
% 
%     convs: [30x7 double]
%       exs: [30x7 double]
%       its: [30x7 double]
%       nbs: [30x7 double]
%      numf: [30x7 double]
%      numg: [30x7 double]
%     times: [30x7 double]
% dataProb = 
% 
%          ms: [409x2 double]
%       msRks: [409x2 double]
%       names: {409x1 cell}
%        nnzs: [409x1 double]
%     probIdx: [100x1 double]

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
                        if ti > 1e2;
                            tfrmt = '%.0f';
                        end
                tblr = [tblr, sprintf(['%i &',tfrmt,'%s'],...
                            round(it),ti,lastc)];  %#ok<AGROW>
%                 tblr = [tblr, sprintf('%i & %.2g %s',...
%                         round(dataEX.its(i,j)),dataEX.times(i,j),lastc)];  %#ok<AGROW>
            end
        
        end
        
%         if j~=skipIdx
%         
%             tblr = [tblr, sprintf('%i & %i & %.2g &',...
%                     dataEX.exs(i,j),round(dataEX.nits(i,j)),dataEX.times(i,j))];  %#ok<AGROW>
%         
%         end
        
    end
    
    table{i} = tblr;%[tblr,'\\'];
    
end

% display(table)