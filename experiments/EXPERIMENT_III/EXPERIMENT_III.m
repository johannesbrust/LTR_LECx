%%---------------------- Experiment III ---------------------------------%%
% This is Experiment 1 in the article
%
% 'Large-Scale Optimization with Linear Equality Constraints'
% J.J. Brust, R.F. Marcia, C.G. Petra, M.A. Saunders (20/21)
%
% Comparison of solvers on the Rosenbrock function, subject to linear
% equality constraints from problems in the SuiteSparse Matrix Collection:
% https://sparse.tamu.edu/
%
% This script is based on a previous version started in 2018.
%--------------------------------------------------------------------------
% Initial version: J.B., 07/03/18
%
%     Change log:
%     07/21/17, J.B., Computes solutions with proposed solver
%     "LMTR_EIG_inf_2_DENSE_G1_CONST_V2" and Matlab's large-scale
%     trust-region optimization solver "fmincon".
%     08/24/17, J.B., Comparison of two versions of the method. [LMTR_..._V2]
%     uses a shape-changing norm on g, whereas [LMTR_..._V3] applies the same
%     norm to an oblique projection of g.
%     08/28/17, J.B., Computations using a modified solver based on a
%     projected gradient stopping condition [LMTR_..._V4].
%     08/28/17, J.B., Computations using a modified solver based on 
%     approximations to the Lagrangian [LMTR_..._V5].
%     09/07/17, J.B., Test on the l-2 based constrained solver
%     10/17/17, J.B, Compact representation and eigenvalues of B^{-1}W solver
%     02/05/18, J.B, Minor edits to prepare numerical results for manuscript.
%     06/05/18, J.B., Use modifed solver that computes performance ratio
%     'rho' to accept steps. Tighter convergence bounds, too. 
%     07/02/18, J.B., This version allows for comparison for large-scale n.
%     07/03/18, J.B., Test on the rosen objective with linear
%     equality constaints selected from netlib.
%     07/04/18, J.B., Comparison with other solvers. RSQP, fmincon-interior,
%     possibly fmincon-trust-region
%     07/19/18, J.B., Preparation for release
%     02/04/19, J.B., Inclusion of Ipopt  
%     10/15/20, J.B., Extension of experiment for new solvers

clc;
clear;

timeEX = tic;

addpath(genpath('../../main'));
addpath(genpath('../../auxiliary'));
addpath(genpath('../../solvers/SPQR/MATLAB'));
addpath(genpath('../../ssget'));

if ismac == 1
    addpath(genpath('../../solvers/IPOPT_UP'));
else
    addpath(genpath('../../solvers/ipopt'));
end

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','..','/data/');
figpath     = fullfile(currentpath,'..','..','/figs/');

rng(090317);

fprintf('---------------- EXPERIMENT III EXT ------------------------\n');

%%------------------ SuiteSparse Matrix Collection Problems --------------%
% Test matrices
%data = load('A_ARWHEAD');
data = load('EVEN_N_PROBS');
probIdx = data.probIdx;

%nprob = length(probIdx);
nprob = 50; % 30

msRks = data.msRks;

numn      = nprob;%length(problems);
numsol    = 7; % 3 SC + 3 L2 + IPOPT

maxM      = 2500; % Maximum size of constraints for solvers not meant for
                  % large problems, i.e., LTRSC_LEC_V1, LTRL2_LEC_V1
                  
%%------------------- Generate data stores ----------------------------- %%

exs             = zeros(numn,numsol);
convs           = zeros(numn,numsol);
nbs             = zeros(numn,numsol);
its             = zeros(numn,numsol);
times           = zeros(numn,numsol);
numf            = zeros(numn,numsol);
numg            = zeros(numn,numsol);

rks             = zeros(numn,1);

%%------------------ Solver parameters --------------------------------- %%
                                 
% %% RSQP 
% options_rsqp            = optimset2('GradObj','on', ...
%                                     'Display','off',...
%                                     'MaxIter', 1e6);

%% L2-Const, SC-Const
options_const.storedat  = 0;
options_const.hastrrad  = 1;
options_const.ctol      = 1e-5;
options_const.btol      = 5e-8;
options_const.dflag     = 0;
options_const.gtol      = 1e-5;
options_const.m         = 5;

%                                     
%% Ipopt
options_ipopt.ipopt.jac_c_constant        = 'yes';
options_ipopt.ipopt.hessian_approximation = 'limited-memory';
options_ipopt.ipopt.mu_strategy           = 'adaptive';
options_ipopt.ipopt.tol                   = 1e-5;
options_ipopt.dual_inf_tol                = 1e-5;
options_ipopt.ipopt.print_level           = 0;
options_ipopt.limited_memory_max_history  = 5;

% SPQR
optsSPQR.Q = 'Householder';

% Convergence tolerance: If nomr(Proj(g)) < ctol and norm(Ax-b) < btol,
% then converged.
ctol                    = 1.e-5;                                    
                                    
%%----------------------- Objective function --------------------------- %%
fun             = @(x)( rosen(x) );   
grad            = @(x)( rosen_ipopt_grad(x));

% 8 spaces reserved for Data display
fprintf('\n**********************\nRunning Solver comparison\n');
fprintf('Row 1: [m,n,rank(A)] \nRow 2: Comput. times per solver \n');
% fprintf(['Every Problem is represented with 2 rows;\n',...
%         'Row 1: [m,n,rank(A)], Row 2: Comput. times per solver. \n']);
fprintf('**********************\n');
fprintf(['m  \t n  \t RNK  \t TR2     \t TR2H    \t',...
            'TR2L    \t TR1     \t TR1H    \t TR1L    \t IPOPT   \n']);

for i = 1:nprob
    % Loading problem and some data    
    Prob= ssget(probIdx(i));
    A   = Prob.A;
    
    [m,n] = size(A);
    %rkA = msRks(i,2);
    
    % Feasible initial point
    x0 = zeros(n,1);
    x0(10) = 5;    
    b = A*x0;
    
    % Computing a factorization for later comparisons and buffer of zeros
    [Q,~,~,infoSP]=spqr(A',optsSPQR);
    rkA = infoSP.rank_A_estimate;
    zrk = zeros(rkA,1);

    if issparse(A); sA = A; else sA = sparse(A); end;
    
    const       = @(x)( const_quad_arg(x,A,b));
    const_ipopt = @(x)( const_quad_arg_ipopt(x,A,b));
    jac         = @(x)( sA );
    
    btol= 1e-7; %norm(A*x0-b)*10;
    
    fprintf(['%i  \t %i  \t %i  \t %s  \t %s  \t %s   \t %s   \t',...
        '%s    \t %s    \t %s    \n'],m,n,rkA,'--------','--------',...
        '--------','--------','--------','--------','--------'); 
    
    %--------------------- Parameters (based on problem data) -----------%%
    
    %% TR1, TR2
    options_const.btol          = btol;
    trradb                      = norm(x0);
    options_const.trradb        = trradb;   
    options_const.whichConv     = 2; % Infinity norm
    options_const_l2            = options_const;
    options_const_l2.maxitroot  = 10;
    options_const_l2.epsroot    = 1e-5;
    
    % Ipopt
    funcs_ipopt.objective           = fun;
    funcs_ipopt.gradient            = grad;
    funcs_ipopt.constraints         = const_ipopt;
    funcs_ipopt.jacobian            = jac;
    funcs_ipopt.jacobianstructure   = jac;

    zm                              = zeros(m,1);
    options_ipopt.cl                = zm;
    options_ipopt.cu                = zm;
    
    sidx                = 0;
    %%-------------------- Solver calls --------------------------------%%
    
    %% TR2
    sidx            = sidx + 1;
    if m < maxM
        [x,~,outinfo]   = LTRSC_LEC_V1(fun,const,x0,options_const); % V6

        [~,g]           = fun(x);

        ztmp    = spqr_qmult(Q,g,0);
        ztmp1   = [zrk;ztmp(rkA+1:end)];
        gp      = spqr_qmult(Q,ztmp1,1);
        
        err             = norm(gp,inf);
        %ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
        %err             = ngp/norm(x);    

        nb              = norm(A*x-b);
    else
        err = NaN;
        nb  = NaN;
        outinfo.numit = NaN;
        outinfo.tcpu = NaN;
        outinfo.numg = NaN;
        outinfo.numf = NaN;
    end
    
    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% TR2H
    sidx            = sidx + 1;
    options_const.whichAsolve = 3;
    
    options_const.dflag = 0;
    
    [x,~,outinfo]   = LTRSC_LEC_V2(fun,const,x0,options_const); % V6
    [~,g]           = fun(x);

    options_const.dflag = 0;
    
    ztmp    = spqr_qmult(Q,g,0);
    ztmp1   = [zrk;ztmp(rkA+1:end)];
    gp      = spqr_qmult(Q,ztmp1,1);
    
    err  = norm(gp,inf);
    
    %err             = norm(g - A'*((A*A')\(A*g)),inf);
        %ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
        %err             = ngp/norm(x);    

    nb              = norm(A*x-b);

    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% TR2L
    sidx            = sidx + 1;
    options_const.whichAsolve = 4;
    [x,~,outinfo]   = LTRSC_LEC_V2(fun,const,x0,options_const); % V6
    [~,g]           = fun(x);
    
    ztmp    = spqr_qmult(Q,g,0);
    ztmp1   = [zrk;ztmp(rkA+1:end)];
    gp      = spqr_qmult(Q,ztmp1,1);
    
    err  = norm(gp,inf);
    
    %err             = norm(g - A'*((A*A')\(A*g)),inf);
        %ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
        %err             = ngp/norm(x);    

    nb              = norm(A*x-b);

    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% TR1
    sidx            = sidx + 1;
    if m < maxM        
        options_const_l2.dflag = 0;
        [x,~,outinfo]   = LTRL2_LEC_V1(fun,const,x0,options_const_l2);    % _V1

        [~,g]           = fun(x);
        
        ztmp    = spqr_qmult(Q,g,0);
        ztmp1   = [zrk;ztmp(rkA+1:end)];
        gp      = spqr_qmult(Q,ztmp1,1);

        err  = norm(gp,inf);
        
    %     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
    %     err             = ngp/norm(x);    
         nb              = norm(A*x-b);
    else
        err = NaN;
        nb  = NaN;
        outinfo.numit = NaN;
        outinfo.tcpu = NaN;
        outinfo.numg = NaN;
        outinfo.numf = NaN;
    end
    
    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;

    %% TR1H
    sidx            = sidx + 1;
    %options_const_l2.whichConv = 2;
    options_const_l2.whichAsolve = 3;
    [x,~,outinfo]   = LTRL2_LEC_V2(fun,const,x0,options_const_l2);    % _V1

    [~,g]           = fun(x);
        
    ztmp    = spqr_qmult(Q,g,0);
    ztmp1   = [zrk;ztmp(rkA+1:end)];
    gp      = spqr_qmult(Q,ztmp1,1);

    nb 	= norm(A*x-b);
    
    err  = norm(gp,inf);

    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;

    %% TR1L
    sidx            = sidx + 1;
    %options_const_l2.whichConv = 2;
    options_const_l2.whichAsolve = 4;
    [x,~,outinfo]   = LTRL2_LEC_V2(fun,const,x0,options_const_l2);    % _V1

    [~,g]           = fun(x);
        
    ztmp    = spqr_qmult(Q,g,0);
    ztmp1   = [zrk;ztmp(rkA+1:end)];
    gp      = spqr_qmult(Q,ztmp1,1);

    err  = norm(gp,inf);
    nb 	= norm(A*x-b);

    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% Ipopt
    sidx                    = sidx + 1;
    %if m < maxM
        [x,info]                = ipopt(x0,funcs_ipopt,options_ipopt);  

        [f,g]                   = fun(x);

        ztmp    = spqr_qmult(Q,g,0);
        ztmp1   = [zrk;ztmp(rkA+1:end)];
        gp      = spqr_qmult(Q,ztmp1,1);

        err  = norm(gp,inf);
        
        nb 	= norm(A*x-b);

        numf(i,sidx)            = -1;
        numg(i,sidx)            = -1;
        
%         ngp                     = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%         err                     = ngp/norm(x);    
%         nb                      = norm(A*x-b);
%     else
%         err = NaN;
%         nb  = NaN;
%         info.numit = NaN;
%         info.cpu = NaN;
%         info.numg = NaN;
%         info.numf = NaN;
%     end
    
    exs(i,sidx)             = (err < ctol) && (nb < btol);
    convs(i,sidx)           = err;
    nbs(i,sidx)             = nb;
    its(i,sidx)             = info.iter;
    times(i,sidx)           = info.cpu;
    
     fprintf(['%i  \t %i  \t %i  \t %.2e \t %.2e \t %.2e \t %.2e \t',...
         '%.2e \t %.2e \t %.2e \n'],m,n,rkA,times(i,1),times(i,2),times(i,3),...
            times(i,4),times(i,5),times(i,6),times(i,7));

end

%%-------------------- Comparison plots ---------------------------------%%
       
leg={   'TR2',...
        'TR2H',...
        'TR2L',...
        'TR1',...       
        'TR1H',...
        'TR1L',...        
        'IPOPT'};
                                            
types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

indAlg          = [1 2 3 4 5 6 7]; % 5 6

perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_III_EXT.eps'));


perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_III_EXT.eps'));

save(fullfile(datapath,'EXPERIMENT_III_EXT'),'exs','convs','nbs','its','times','numf','numg');

close ALL;

% Restore warning settings
warning(wtest);

timeEX = toc(timeEX);
