%-------------------- table_EX -------------------------------------------%
%
% Script to store data from Experiment into a table for the
% manuscript
%
% The table has 61 rows and 2 + 2 * 7 = 17 columns (b/c one small problem
% is skipped)
%
%-------------------------------------------------------------------------%
% 10/16/20, J.B.
% 10/21/20, J.B., problem mask
% 10/27/20, J.B., Table for second experiment
% 11/02/20, J.B., Modify time result format depending on size

% Load data from experiment
datapath = fullfile(pwd,'..','..','/data/EXPERIMENT_IV_EXT');
dataEX = load(datapath);

probpath = fullfile(pwd,'..','/..','/auxiliary/');

% Load problem names
fid     = fopen(fullfile(probpath,'cutest_list.txt'),'r');
%dataProb = load('EVEN_N_PROBS');
%probMask = dataProb.probMask;

nrow    = 61; % 30
%ncol    = 4 + 3*9;
table = cell(nrow,1);
nsol    = 7;
skipIdx = 0; % Don't include reprojected solver

% Data stored in
% dataEX = 
% 
% 'ex','numit','t_aver','numf','numg','params','tract', 'nms'

i = 0;
tidx = 0;
tline = fgets(fid);
while ischar(tline)
    tline = fgets(fid);
    
    if (~strcmp(tline(1),'%'))  && (ischar(tline))
        
        i = i + 1;
        
        n = dataEX.nms(i,2);
        m = dataEX.nms(i,1);
        
        if 500 < n
            
            tidx = tidx + 1;
            
            % Process data
            tline_ = tline;
            pname = strtrim(tline);
%             pname_L = dataProb.names{probMask(i)};
%             pname_SP= strsplit(pname_L,'/');
%             pname = pname_SP{2};
            
%             pname = strrep(pname,'_','\_');
%             m = dataProb.ms(probMask(i),1);
%             n = dataProb.ms(probMask(i),2);
%             rnk = dataProb.msRks(probMask(i),2);
%             dens = (dataProb.nnzs(probMask(i))/(m*n));
            
            %tblr = sprintf('$\\texttt{%s}$ & %i/%i & %i/%.1g &',pname,m,n,rnk,dens);
            tblr = sprintf('$\\texttt{%s}$ & %i/%i &',pname,m,n);
            %dataEX.conds(i));
            for j = 1:nsol
                
                lastc = '&';
                if j == nsol
                    lastc ='\\';
                end
                
                it = dataEX.numit(i,j);
                ti = dataEX.t_aver(i,j);
                
                if j~=skipIdx
                    % Different row depending on outcome
                    ex = dataEX.ex(i,j); % Either 1,-1,-2
                    % Nonconvergence condition
                    cond1 = ((ex < 1) && (j < nsol) && (it < 100001))||...
                            ((ex < 1) && (j == nsol) && (it < 3000));
                    % Max iter condition
                    cond2 = ((ex < 1) && (j < nsol) && (it == 100001))||...
                            ((ex < 1) && (j == nsol) && (it == 3000));    
                    if cond1 == 1
                        tblr = [tblr, sprintf('$\\texttt{NC}^{\\dagger}$ & $\\texttt{NC}$ %s',lastc)]; %#ok<AGROW>
                        %tblr = [tblr, sprintf('$\\texttt{N/A}^{*}$ & $\\texttt{N/A}$ %s',lastc)]; %#ok<AGROW>
                    elseif cond2 == 1
                        tblr = [tblr, sprintf('$\\texttt{MX}^{\\dagger}$ & $\\texttt{MX}$ %s',lastc)]; %#ok<AGROW>
                    else
                        tfrmt = '%.2g';
                        if ti > 1e2;
                            tfrmt = '%.0f';
                        end
                        tblr = [tblr, sprintf(['%i &',tfrmt,'%s'],...
                            round(it),ti,lastc)];  %#ok<AGROW>
                    end
                    
                end
                
                %         if j~=skipIdx
                %
                %             tblr = [tblr, sprintf('%i & %i & %.2g &',...
                %                     dataEX.exs(i,j),round(dataEX.nits(i,j)),dataEX.times(i,j))];  %#ok<AGROW>
                %
                %         end
                
            end
            
            table{tidx} = tblr;%[tblr,'\\'];
            
        end
        
    end
end

fclose(fid);
% display(table)