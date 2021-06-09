%-------------------- Experiment_V_EXT -----------------------------------%
% From the article 'Large-Scale Optimization With Linear Equality Constraints
% using Reduced Compact Representation' J.J. Brust, R.F. Marcia,
% C.G. Petra and M.A. Saunders (2021)
%
% Test problems on large CUTEst problems, with linear equality constraints.
% Problems that have intrinsically bounds, are treated as if there are
% unbounded. Problems 19-31 are shifted by rho.*I (with rho=10)
%
% Comparison between two proposed solvers and three external solvers
% (depending on context)
%
%--------------------------------------------------------------------------
% Initial version: 09/14/17, J.B.,

% 06/05/18, J.B., Use modifed solver that computes performance ratio
%    'rho' to accept steps. Tighter convergence bounds, too.
% 06/27/16, J.B., Comparisons of proposed solvers with fminunc and RSQP
% 07/05/18, J.B., Comparison between the two proposed solvers, and fmincon
% interior-point solver.
% 07/08/18, J.B., Removal of fmincon solver
% 07/19/18, J.B., Preparation for release
% 02/04/19, J.B., Inclusion of Ipopt
% 10/14/20, J.B., Initial setup for test of extended methods
% 10/22/20, J.B., Test of comparing 7 solvers
% 05/25/21, J.B., Additional test problems
% 05/27/21, J.B., New test problems, preparation of the experiment
% 05/28/21, J.B., Update of experiment
% 06/02/21, J.B., Experiment with selected new problems
% 06/03/21, J.B., Limiting the size of the number of equality constraints
%                   for TR1, TR2
% 06/07/21, J.B., Update of experiment (skip problem BLOWEYC for TR1,TR2)
% 06/08/21, J.B., Take out all BLOWEY* problems
% 06/09/21, J.B., Preparation for release
%--------------------------------------------------------------------------
%NOTE(S): This Experiment requires CUTEst to be installed

clc;
clear;

addpath(genpath('../../main'));
addpath(genpath('../../auxiliary'));
addpath(genpath('../../netlib/readmps'));
addpath(genpath('../../solvers/SPQR/MATLAB'));

if ismac == 1
    addpath(genpath('../../solvers/IPOPT_UP'));
else
    addpath(genpath('../../solvers/ipopt'));
end

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','/..','/data/');
figpath     = fullfile(currentpath,'..','/..','/figs/');
probpath    = fullfile(currentpath,'..','/..','/auxiliary/');
%probpath    = fullfile(currentpath);

rng(090317);

fprintf('---------------- EXPERIMENT Va ------------------------\n');
tEX = tic;

%%----------------------- Parameters ----------------------------------%

% Set input parameters:
% %% Proposed solvers
% params=struct;
% params.m = 5;  % number of L-BFGS updates
% params.gtol = 1e-5;  % exit if ||g||_2<gtol*max(1,||x||_2)
% params.ranktol = 1e-7;  % tolerance for establishing rank of V
% params.dflag = 0;  % display parameter, 1 if to display information per iteration
% params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
% params.ftol=1e-11;  % tolerance on relative function reduction
%
% params.storedat =0;
% params.hastrrad =1;
% params.ctol     = 5e-5;
% params.btol     = 1e-10;
% params.dflag    = 0;
% params.maxit    = 100000;

%% L2-Const, SC-Const
options_const.storedat  = 0;
options_const.hastrrad  = 1;
options_const.ctol      = 1e-5;
options_const.btol      = 5e-8;
options_const.dflag     = 0;
options_const.gtol      = 1e-5;
options_const.m         = 5;
options_const.whichInit = 1;
options_const.whichConv = 1;
%options_const.maxit     = 50;

%% fmincon
% The extended test does not include fmincon
% interior-point algorithm
% options_fmin_int        = optimoptions('fmincon',...
%                                         'GradObj',      'on',...
%                                         'Algorithm',    'interior-point', ...
%                                         'Hessian',      'lbfgs',...
%                                         'MaxIter',      1e5, ...
%                                         'MaxFunEvals',  1e6, ...
%                                         'TolX',         1e-10,...
%                                         'SubproblemAlgorithm', 'ldl-factorization'); % 'SubproblemAlgorithm', 'ldl-factorization'

%% Ipopt
options_ipopt.ipopt.jac_c_constant        = 'yes';
options_ipopt.ipopt.hessian_approximation = 'limited-memory';
options_ipopt.ipopt.mu_strategy           = 'adaptive';
options_ipopt.ipopt.tol                   = 1e-5;
%options_ipopt.dual_inf_tol                = 1e-5;
options_ipopt.ipopt.print_level           = 0; % 5
%options_ipopt.limited_memory_max_history  = 5;
options_ipopt.ipopt.max_iter              = 100000; % 5

CUTEst_init  % initialize CUTEr, see appropriate documentation
%fid     = fopen(fullfile(probpath,'probs_EX_V.txt'),'r');
fid     = fopen(fullfile(probpath,'cutest_list_EX_V.txt'),'r');

sfil    = 'TestResults';

% Initialize storage for output information
numRuns         = 1; % 3
numAlgorithms   = 7; % 7
numProblems     = 31;
ex              = zeros(numProblems,numAlgorithms);
numf            = zeros(numProblems,numAlgorithms);
numg            = zeros(numProblems,numAlgorithms);
numit           = zeros(numProblems,numAlgorithms);
tcpu            = zeros(numProblems,numRuns,numAlgorithms);
t_aver          = zeros(numProblems,numAlgorithms);
tract           = zeros(numProblems,numAlgorithms);
numrst          = zeros(numProblems,numAlgorithms);

nms             = zeros(numProblems,2);
numDenseA       = zeros(numProblems,1); % Storing rm dense cols.

p=1;

pmax = numProblems; % 1/numProblems; % testing/All problems
maxM = 2500; % Maximum size of constraints for TR1 and TR2

% Options SPQR
optsSPQR.Q = 'Householder';

tline = fgets(fid);

s = 0; % Solver index
pIdx = 0; % Problem index from file

% Global shift parameter for problems selected problems
lpIdx = 19;
upIdx = 31;
global rho;
rho = 10;

% Problem(s) to skip for TR1, TR2 because they take very long.
% (These methods were not developed for large-scale problems)
clist ={'BLOWEYA','BLOWEYB',...
        'BLOWEYC'};

while ischar(tline)
    
    pIdx = pIdx + 1;
    
    if (~strcmp(tline(1),'%'))  && (ischar(tline)) && p <= pmax
        
        if isunix == 1 && ismac == 0
            eval(['!runcutest -p matlab -D ' tline]);
        else
            tline_ = tline;
            tline = strtrim(tline);
            cmdcutest   = ['cutest2matlab_osx $MASTSIF/' tline];
            unix(cmdcutest);
        end
        
        prob            = cutest_setup();
        x0              = prob.x;
        params.trradb   = max(norm(x0),1);
        n               = size(x0,1);
        
        % CUTEST constraint
        [c,Ad] = cutest_cons(x0);
        A = sparse(Ad);
        mm = size(A,1);
        b = -cutest_cons(zeros(n,1));
        
        nms(p,1) = n;
        nms(p,2) = mm;
        
        % Making problem feasible if neccessary
        if norm(c) > 1e-10
            %x0 = x0 - spqr_solve(A,c);
            x0 = x0 - spqr_solve(A,c,struct('solution','min2norm'));
        end
        
        cons = @(x)( const_quad_arg(x,A,b));
        
        % Shifts the problem Hessian by "rho"
        if (lpIdx <= pIdx) && (pIdx <= upIdx)
            obj = @cutest_fun_shifted;
            grad = @(x)(cutest_fun_shifted(x,'G'));
        else
            obj = @cutest_fun;
            grad = @(x)(cutest_fun(x,'G'));
        end
        
        sA              = sparse(A);
        const_ipopt     = @(x)( const_quad_arg_ipopt(x,A,b));
        jac             = @(x)( sA );
        
        % Ipopt
        funcs_ipopt.objective           = obj;
        funcs_ipopt.gradient            = grad;
        funcs_ipopt.constraints         = const_ipopt;
        funcs_ipopt.jacobian            = jac;
        funcs_ipopt.jacobianstructure   = jac;
        
        zm                              = zeros(mm,1);
        options_ipopt.cl                = zm;
        options_ipopt.cu                = zm;
        
        options_const_l2            = options_const;
        options_const_l2.maxitroot  = 10;
        options_const_l2.epsroot    = 1e-5;
        
        % TR1
        s=1;
        if mm < maxM && sum(strcmp(tline,clist))==0
            [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
                runAlgorithm(@LTRSC_LEC_V1,obj,cons,x0,options_const,numRuns); % LTRSC_LEC_V1
        else
            ex(p,s) = NaN;
            numf(p,s) = NaN;
            numg(p,s) = NaN;
            numit(p,s) = NaN;
            tcpu(s,:,p) = NaN;
            tract(p,s) = NaN;
        end
        
        % TR1x1
        s=s+1;
        options_const.whichAsolve = 3;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LTRSC_LEC_V2,obj,cons,x0,options_const,numRuns); % LTRL2_LEC_V1
        
        % TR1x2
        s=s+1;
        options_const.whichAsolve = 4;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LTRSC_LEC_V2,obj,cons,x0,options_const,numRuns); % LTRL2_LEC_V1
        
        % TR2
        s=s+1;
        options_const_l2.whichAsolve = 1;
        if mm < maxM && sum(strcmp(tline,clist))==0
            [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
                runAlgorithm(@LTRL2_LEC_V1,obj,cons,x0,options_const_l2,numRuns); % LTRL2_LEC_V1
        else
            ex(p,s) = NaN;
            numf(p,s) = NaN;
            numg(p,s) = NaN;
            numit(p,s) = NaN;
            tcpu(s,:,p) = NaN;
            tract(p,s) = NaN;
        end
        
        % TR2x1
        s=s+1;
        options_const_l2.whichAsolve = 3;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LTRL2_LEC_V2,obj,cons,x0,options_const_l2,numRuns); % LTRL2_LEC_V1
        % TR2x2
        s=s+1;
        options_const_l2.whichAsolve = 4;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LTRL2_LEC_V2,obj,cons,x0,options_const_l2,numRuns); % LTRL2_LEC_V1
        %
        %         % fmincon
        % %        s = s+1;
        % %         for ir = 1:numRuns
        % %            tStart = tic;
        % %            options_fmin_int.OutputFcn       = @(x,oV,state)(outputFcn_fmincon(x,oV,state,tStart));
        % %
        % %
        % %            [x_fmin,f_fmin,ex_fmin,out_fmin] = fmincon(@cutest_fun,x0,[],[],A,b,...,
        % %                                             [],[],[],options_fmin_int);
        % %
        % %            time_alg         = toc(tStart);
        % %            tcpu(s,ir,p)     = time_alg;
        % %            ex(p,s)          = ex_fmin;
        % %            numf(p,s)        = out_fmin.funcCount;
        % %            numg(p,s)        = out_fmin.funcCount;
        % %            numit(p,s)       = out_fmin.iterations;
        % %         end
        %
        
        % Ipopt
        s = s+1;
        % exit : 0,1,2 acceptable
        for ir = 1:numRuns
            tStarti = tic;
            %          funcs_ipopt.iterfunc       = @(x,oV,state)(...
            %                                            outputFcn_ipopt(x,oV,state,tStarti));
            
            [x,info]                = ipopt(x0,funcs_ipopt,options_ipopt);
            
            if ( (-1 < info.status) && ( info.status < 3 )), ex_ipopt = 1; else ex_ipopt = -1;end;
            
            time_alg         = toc(tStarti);
            tcpu(s,ir,p)     = info.cpu;
            ex(p,s)          = ex_ipopt;
            numf(p,s)        = -1;
            numg(p,s)        = -1;
            numit(p,s)       = info.iter;
        end
        
        % Average CPU time
        if p==1 && numRuns > 2
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        elseif numRuns == 2
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        else
            for si=1:s
                t_aver(p,si) = tcpu(si,1:1,p);
            end
        end
        
        cutest_terminate();
        
        delete( '*.d',...
            '*.o',...
            '*.dylib',...
            '*.f',...
            'mcutest.*');
        
        p=p+1;
    end
    
    tline = fgets(fid);
    
end

leg={   'TR1',...
    'TR1x1',...
    'TR1x2',...
    'TR2',...
    'TR2x1',...
    'TR2x2',...
    'IPOPT'};

types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

indAlg          = [1 2 3 4 5 6 7];


perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_V_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_V_EXT.eps'));

% Comparison between proposed solvers only (TR1)
indAlg          = [1 2 3];

perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_SEL_TR1_EX_V_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_SEL_TR1_EX_V_EXT.eps'));

% Comparison between proposed solvers only (TR2)
indAlg          = [4 5 6];

perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_SEL_TR2_EX_V_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_SEL_TR2_EX_V_EXT.eps'));

save(fullfile(datapath,'EXPERIMENT_V_EXT'),'ex','numit','t_aver','numf',...
    'numg','params','tract','nms','numDenseA');

close ALL;

% Restore warning settings
warning(wtest);

tEX = toc(tEX); % Time to run experiment