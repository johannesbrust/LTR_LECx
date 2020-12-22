%%---------------------- Experiment III ---------------------------------%%
%%(Check for IPOPT behavior on Experiment 1)
% From the article 'Large-Scale Quasi-Newton Trust-Region Methods
% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
% C.G. Petra
%
% Comparison of solvers on the Rosenbrock function, subject to linear
% equality constraints from netlib problems: http://www.netlib.org/lp/data/
%
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
%     11/06/20, J.B., Monitoring IPOPT (problem 1 "beacxc")

clc;
clear;

timeEX = tic;

addpath(genpath('../../main'));
addpath(genpath('../../auxiliary'));
%addpath ../solvers/rsqp
%addpath ../netlib/readmps
addpath(genpath('../../solvers/IPOPT_UP'));
addpath(genpath('../../solvers/SPQR/MATLAB'));
addpath(genpath('../../ssget'));

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','..','/data/');
figpath     = fullfile(currentpath,'..','..','/figs/');
%probpath    = fullfile(currentpath,'..','..','/netlib/');

rng(090317);

fprintf('---------------- EXPERIMENT III EXT ------------------------\n');

%%------------------ SuiteSparse Matrix Collection Problems --------------%
% Test matrices
%data = load('A_ARWHEAD');
data = load('EVEN_N_PROBS');
probIdx = data.probIdx;

%nprob = length(probIdx);
nprob = 50;%50; % 30; % 1

msRks = data.msRks;
% problems = {'fit1d.mps',...
%             'fit2d.mps',...
%             'd6cube.mps',...
%             'scsd1.mps',...
%             'scsd6.mps',...
%             'scsd8.mps'};

numn      = nprob;%length(problems);
numsol    = 7; % 3 SC + 3 L2 + IPOPT

maxM      = 2500; % Maximum size of constraints
%%------------------- Generate data stores ---------------------------- %%

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

%% fmincon
% interior-point algorithm
% options_fmin_int        = optimoptions('fmincon',...   
%                                         'GradObj',      'on',...
%                                         'Display',      'off',...
%                                         'Algorithm',    'interior-point', ...
%                                         'Hessian',      {'lbfgs',5},...
%                                         'MaxIter',      1e6, ...
%                                         'MaxFunEvals',  1e6, ...
%                                         'TolX',         1e-10);
% % trust-region-reflective algorithm
% options_fmin_tr        = optimoptions('fmincon',...   
%                                         'GradObj',      'on',...
%                                         'Display',      'off',...
%                                         'Algorithm',    'trust-region-reflective', ...                                        
%                                         'MaxIter',      1e6, ...
%                                         'MaxFunEvals',  1e6);
%                                     
%% Ipopt
options_ipopt.ipopt.jac_c_constant        = 'yes';
options_ipopt.ipopt.hessian_approximation = 'limited-memory';
options_ipopt.ipopt.mu_strategy           = 'adaptive';
options_ipopt.ipopt.tol                   = 1e-5;
options_ipopt.dual_inf_tol                = 1e-5;
options_ipopt.ipopt.print_level           = 5;%0;
options_ipopt.limited_memory_max_history  = 5;

% SPQR
optsSPQR.Q = 'Householder';

% Convergence tolerance: If nomr(Proj(g)) < ctol and norm(Ax-b) < btol,
% then converged.
ctol                    = 1.e-5;                                    
                                    
%%----------------------- Objective function --------------------------- %%
fun             = @(x)( rosen(x) );   
grad            = @(x)( rosen_ipopt_grad(x));

fprintf('\n**********************\nRunning Solver comparison\n');
%fprintf('n\t time RSQP\t time SC\t time L2\t time  fmincon-I\t fmincon-TR\n');
%fprintf('n\t m\t rnk\t t-TRSC\t t-TRL2\t t-fmin \t t-Ipopt\n');
fprintf(['m  \t n  \t RNK  \t TR1    \t TR1x1  \t',...
            'TR1x2  \t TR2  \t TR2x1  \t TR2x2 \t IPOPT  \n']);
% fprintf(['m  \t n  \t RNK  \t TR1    \t TR1x1  \t',...
%             'TR2  \t TR2x1  \t IPOPT  \n']);       
%nprob = 1;
clist={'beacxc';...
    'lp_cre_d';...
    'lp_fit2d';...
    'lp_scsd1'};
%clist={'lp_scsd1'};
    
for i = 1:nprob
    % Loading problem and some data    
    Prob= ssget(probIdx(i));
    A   = Prob.A;
    
    [m,n] = size(A);
    %rkA = msRks(i,2);
    
%     name    = problems{i};
%     prob    = readmps(fullfile(probpath,name));
%     
%     A       = prob.A;
%     [mn,n]  = size(A);
%     
%     rnk     = rank(full(A));                      
%     
%     b0      = randn(n,1);
%     b       = A*b0;
%     x0      = A'*((A*A')\b);

    % Process problem only if not in list
    pname_L = Prob.name;
    pname_SP= strsplit(pname_L,'/');
    pname = pname_SP{2};
    
    if sum(strcmp(pname,clist))==0
        continue; % Skip problem if not in clist
    end
    
    
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
    
%     fprintf(['%i  \t %i  \t %i  \t %s  \t %s  \t %s   \t %s   \t',...
%         '%s    \n'],m,n,rkA,'------','------','------','------',...
%         '------'); 
    
    fprintf(['%i  \t %i  \t %i  \t %s  \t %s  \t %s   \t %s   \t',...
        '%s    \t %s    \t %s    \n'],m,n,rkA,'------','------','------','------',...
        '------','------','------'); 
    
    %--------------------- Parameters (based on problem data) -----------%%
    
    %% L2-Const, SC-Const    
    options_const.btol          = btol;
    trradb                      = norm(x0);
    options_const.trradb        = trradb;   
    options_const.whichConv     = 2; % Infinity norm
    options_const_l2            = options_const;
    options_const_l2.maxitroot  = 10;
    options_const_l2.epsroot    = 1e-5;
    
%     %% fmincon
%     options_fmin_int    = optimoptions(options_fmin_int,...
%                                         'TolCon', btol);
%     options_fmin_tr     = optimoptions(options_fmin_tr,...
%                                         'TolCon', btol);
                                    
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
%     %% RSQP
%     sidx                = sidx + 1;
%     
%     tic;
%     [x,f_rsqp,ex_rsqp,out_rsqp] = ...
%     rsqp(fun,x0,A,b,A,b,[],[],[],options_rsqp);
%     time_rsqp = toc;
%     
%     [f,g]           = fun(x);
%     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%     err             = ngp/norm(x);    
%     nb              = norm(A*x-b);
%     
%     exs(i,sidx)     = (err < 1e-5) && (nb < btol);
%     convs(i,sidx)   = err;
%     nbs(i,sidx)     = nb;
%     its(i,sidx)     = out_rsqp.iterations;
%     times(i,sidx)   = time_rsqp;
%     numf(i,sidx)    = out_rsqp.funcCount;
%     numg(i,sidx)    = NaN;
    
    %% TR1
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
    
    %% TR1x1
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
    
    %% TR1x2
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
    
    %% TR2
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

    %% TR2x1
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

    %% TR2x2
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
    
    %% fmincon-I
%     sidx            = sidx + 1;
%     tic;
%     [x,~,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
%         [],[],[],options_fmin_int);
%     time_fmincon = toc;
%     
%     [f,g]           = fun(x);
%     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%     err             = ngp/norm(x);    
%     nb              = norm(A*x-b);
%     
%     exs(i,sidx)     = (err < ctol) && (nb < btol);
%     convs(i,sidx)   = err;
%     nbs(i,sidx)     = nb;
%     its(i,sidx)     = out_fmin.iterations;
%     times(i,sidx)   = time_fmincon;
%     numf(i,sidx)    = out_fmin.funcCount;
%     numg(i,sidx)    = out_fmin.funcCount;
%     
    %fprintf('%i\t %i\t %i\t %6.4f\t %6.4f\t %6.4f\n',n,mn,rnk,times(i,1),times(i,2),times(i,3));
    
%     %% fmincon-TR
%     sidx            = sidx + 1;
%     tic;
%     [x,f,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
%         [],[],[],options_fmin_tr);
%     time_fmincon = toc;
%     
%     [f,g]           = fun(x);
%     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%     err             = ngp/norm(x);    
%     nb              = norm(A*x-b);
%     
%     exs(i,sidx)     = (err < 1e-5) && (nb < btol);
%     convs(i,sidx)   = err;
%     nbs(i,sidx)     = nb;
%     its(i,sidx)     = out_fmin.iterations;
%     times(i,sidx)   = time_fmincon;
%     numf(i,sidx)    = out_fmin.funcCount;
%     numg(i,sidx)    = out_fmin.funcCount;
%     
%     %%-------------- Display times -------------------------------------%%
%     fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e \t %.3e\n',n,...
%         times(i,1),...
%         times(i,2),...
%         times(i,3),...
%         times(i,4),...
%         times(i,5));    

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
        
%      fprintf(['%i  \t %i  \t %i   \t %.2e \t %.2e \t %.2e \t %.2e \t %.2e ',...
%          '\n'],m,n,rkA,times(i,1),times(i,2),times(i,3),...
%             times(i,4),times(i,5));   

end

%%------------------- Comparison plots ---------------------------------%%

%   
% leg={   'RSQP',...
%         'TR-$(\mathbf{P},\infty)$',...
%         'TR-$\ell_2$',...       
%         'fmincon-I-ldl',...
%         'fmincon-I-cg'};
       
leg={   'TR1',...
        'TR1x1',...
        'TR1x2',...
        'TR2',...       
        'TR2x1',...
        'TR2x2',...        
        'IPOPT'};
% leg={   'TR1',...
%         'TR1x1',...        
%         'TR2',...       
%         'TR2x1',...        
%         'IPOPT'};    
                            
%types.colors    = ['k' 'b' 'r' 'm' 'g']; %'k' 'y'
%types.lines     = {'-.', '-', '-.', '-','.-' }; %'-',   '-'
%types.markers   = ['o' 'o' 'o' 'o' 'o']; %'s' 's'
                
types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

%indAlg          = [1 2 3 4 5]; %

indAlg          = [1 2 3 4 5 6 7]; % 5 6

perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_III_EXT_IPOPTCHK.eps'));


perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_III_EXT_IPOPTCHK.eps'));

save(fullfile(datapath,'EXPERIMENT_III_EXT_IPOPTCHK'),'exs','convs','nbs','its','times','numf','numg');

close ALL;

% Restore warning settings
warning(wtest);

timeEX = toc(timeEX);
