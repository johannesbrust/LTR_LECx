%-------------------- table_EX_IV ----------------------------------------%
%
% Script to store data from Experiment into a table for the
% manuscript
%
% The table has 61 rows and 2 + 2 * 7 = 16 columns (b/c one small problem
% is skipped)
%
%-------------------------------------------------------------------------%
% 10/16/20, J.B.
% 10/21/20, J.B., problem mask
% 10/27/20, J.B., Table for second experiment
% 11/02/20, J.B., Modify time result format depending on size
% 06/07/21, J.B., Patch data of IPOPT with "re-runs"
% 06/09/21, J.B., Preparation for release

% Load data from experiment
datapath = fullfile(pwd,'..','..','/data/');
dataEX = load([datapath,'EXPERIMENT_IV_EXT']);
dataS  = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTED']);
dataSa = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTEDa']);
dataSb = load([datapath,'EXPERIMENT_IV_EXT_IPOPT_RERUN_SELECTEDb']);

% Backup for testing
%dataEXBU = dataEX;

% Patching column of IPOPT
colIP = 7;
for p=1:(size(dataEX.ex,1))
    
    if dataS.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataS.ex(p,colIP);
       dataEX.numit(p,colIP) = dataS.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataS.t_aver(p,colIP);       
       
    elseif dataSa.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataSa.ex(p,colIP);
       dataEX.numit(p,colIP) = dataSa.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataSa.t_aver(p,colIP);
       
    elseif dataSb.ex(p,colIP) == 1
        
       dataEX.ex(p,colIP) = dataSb.ex(p,colIP);
       dataEX.numit(p,colIP) = dataSb.numit(p,colIP);
       dataEX.t_aver(p,colIP) = dataSb.t_aver(p,colIP); 
       
    end
    
end

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
            
            % Sort outcomes
            sit = sort(dataEX.numit(i,:),'ascend');
            sti = sort(dataEX.t_aver(i,:),'ascend');
            
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
                            ((ex < 1) && (j == nsol) && (it < 100001)); % 3000
                    % Max iter condition
                    cond2 = ((ex < 1) && (j < nsol) && (it == 100001))||...
                            ((ex < 1) && (j == nsol) && (it == 100001)); % 3000   
                    if cond1 == 1
                        tblr = [tblr, sprintf('$\\texttt{NC}^{\\dagger}$ & $\\texttt{NC}$ %s',lastc)]; %#ok<AGROW>
                        %tblr = [tblr, sprintf('$\\texttt{N/A}^{*}$ & $\\texttt{N/A}$ %s',lastc)]; %#ok<AGROW>
                    elseif cond2 == 1
                        tblr = [tblr, sprintf('$\\texttt{MX}^{\\dagger}$ & $\\texttt{MX}$ %s',lastc)]; %#ok<AGROW>
                    else
                        tfrmt = '%.2g';
                        
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
                        
%                         if ti > 1e2;
%                             tfrmt = '%.0f';
%                         end
                        
                        
                        tblr = [tblr, sprintf(['%i &',tfrmt,'%s'],...
                            round(it),ti,lastc)];  %#ok<AGROW>
                    end
                    
                end
                               
            end
            
            table{tidx} = tblr;%[tblr,'\\'];
            
        end
        
    end
end

fclose(fid);

% Write to file
fid_w   = fopen([datapath,'table_EX4.txt'],'w');
fprintf(fid_w,'%s \n',table{:});
fclose(fid_w);

% display(table)