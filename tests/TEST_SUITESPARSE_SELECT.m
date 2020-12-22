%--------------- TEST_SUITESPARSE_SELECT ---------------------------------%
%
% Script to select test problems from the SuiteSparse collection of 
% test matrices. 
% Script selects matrices A corresponding to a functon with m = 5000
% variables (i.e., the ARWHEAD functions from CUTEst)
%-------------------------------------------------------------------------%
% 10/15/20, J.B., Modification of a previous version
% 10/21/20, J.B., Adding a mask to retrieve m/n problem information

clc;
clear;

addpath('../ssget');

index = ssget;


% Select rectangular matrices with sizes nl <= n <= nu and m fixed.

mfix    = 5000;
nl      = 500;
nu      = 4000;

% condit =  (index.nrows == mfix) & ...
%           (index.ncols <= nu) & ...
%           (nl <= index.ncols);

% condit =  (index.nrows == mfix) & ...
%           (index.ncols ~= index.nrows);

condit =  (index.ncols > index.nrows) & ...
            (index.ncols > nl);
          
            
% Linear indices    
ids = find(condit);
nids = length(ids);

% Storing problem data
nprob   = nids; %50
msRks   = zeros(nprob,2);
conds   = zeros(nprob,1);
names   = cell(nprob,1);
iids     = zeros(nprob,1);
nnzs    = zeros(nprob,1);
ms      = zeros(nprob,2);
issym   = zeros(nprob,1);
%pidxs   = zeros(nprob,1);

colprev=0;
pidx = 0;
nsel = 100;
probIdx = zeros(nsel,1);
probMask = zeros(nsel,1);

%symtol = 1e-9;

fprintf('----------- Running problem probing ----------- \n');
mess = 'm \t n \t cond(A) \t rank \t nnz \t issym \t prb.time \n';
% for i = 1:nsol    
%     if i==nsol; est = '\n'; else est = '   \t';  end;
%     mess = [mess,'Sol.',num2str(i),est]; %#ok<AGROW>
% end
fprintf(mess);


for i = 1:nprob
    
    pts = tic;
    Prob= ssget(ids(i));
    A   = Prob.A;
    
    [rw,col]   = size(A);
    rk  = sprank(A);
    %rk  = rank(full(A));
    msRks(i,:) = [col,rk];
%     condA = cond(full(A));
%     conds(i) = condA;
    iids(i)   = Prob.id;
    names{i} = Prob.name;
    nnzs(i) = nnz(A);
    ms(i,:) = [rw,col];    
    pte = toc(pts);
    
    %symerr = norm(A-A','fro');
    
    %issym(i) = symerr < symtol;
    
    %messb = '%i \t --- \t --- \t %i \t %i \t %.2e \n';
    messb = '%i \t %i \t ------- \t %i \t %i \t %i \t %.2e \n';
    
    fprintf(messb,rw,col,rk,nnzs(i),issym(i),pte);
    
%     if abs(colprev-col)~=0 %&& issym(i)~=true
%         pidx = pidx+1;
%         probIdx(pidx) = ids(i);
%         if pidx == nsel
%             break;
%         end
%     end
%     colprev = col;
    
    if mod(col,2) == 0 && pte < 1e-1
        pidx = pidx+1;
        probIdx(pidx) = ids(i);
        probMask(pidx) = i;
        if pidx == nsel
          break;
        end
    end
end
% msRks   = zeros(nprob,2);
% conds   = zeros(nprob,1);
% names   = cell(nprob,1);
% iids     = zeros(nprob,1);
% nnzs    = zeros(nprob,1);
% ms      = zeros(nprob,1);
% issym   = zeros(nprob,1);
save('EVEN_N_PROBS','probIdx','msRks','names','nnzs','ms', 'probMask');




%nsym = sum(issym);