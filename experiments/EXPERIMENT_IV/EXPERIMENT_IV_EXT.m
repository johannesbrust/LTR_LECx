%-------------------- Experiment_IV_EXT ---------------------------------%%
% Originally from the article 'Large-Scale Quasi-Newton Trust-Region Methods
% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
% C.G. Petra
% Test problems generated by using a list of unconstrained
% CUTEst problems, and adding linear equality constraints to them. The
% linear constraints are synthetically generated.
%
% Comparison between two proposed solvers and two external solvers
%
% Extended experiment for the manuscript:
% 'Efficient Large-Scale Quasi-Newton Trust-Region Methods With 
% Linear Equality Constraints', 2020
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
% 05/27/21, J.B., Re-run of experiment without iteration function for IPOPT
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

rng(090317);

fprintf('---------------- EXPERIMENT IV_EXT ------------------------\n');
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
options_ipopt.dual_inf_tol                = 1e-5;
options_ipopt.ipopt.print_level           = 0;
options_ipopt.limited_memory_max_history  = 5;                                 
                                    
CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid     = fopen(fullfile(probpath,'cutest_list.txt'),'r');  % read file that contains CUTEr test problems
sfil    = 'TestResults';

% Initialize storage for output information
numRuns         = 1; % 3
numAlgorithms   = 7;
numProblems     = 62;
%numProblems     = 5;
ex              = zeros(numProblems,numAlgorithms);
numf            = zeros(numProblems,numAlgorithms);
numg            = zeros(numProblems,numAlgorithms);
numit           = zeros(numProblems,numAlgorithms);
tcpu            = zeros(numProblems,numRuns,numAlgorithms);
t_aver          = zeros(numProblems,numAlgorithms);
tract           = zeros(numProblems,numAlgorithms);
numrst          = zeros(numProblems,numAlgorithms);

nms             = zeros(numProblems,2);

%mm              = 10; % Number of equality constraints

p=1;

pmax = 62; % All problems

tline = fgets(fid);
while ischar(tline)     
    tline = fgets(fid);       
    
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
        
        %mm = 10;
        mm              = ceil(0.25*n);
        nms(p,1)        = mm;
        nms(p,2)        = n;
        density = 0.1;
        A               = (sprandn(mm,n,density))./params.trradb;
        
        b0              = randn(n,1);    
        b               = A*b0;
        
        x0              = A'*((A*A')\b);
        
        cons            = @(x)( const_quad_arg(x,A,b));
        obj             = @cutest_fun;
        grad            = @cutest_gradient;
        
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
         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
                runAlgorithm(@LTRSC_LEC_V1,obj,cons,x0,options_const,numRuns); % LTRSC_LEC_V1
            
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
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
                runAlgorithm(@LTRL2_LEC_V1,obj,cons,x0,options_const_l2,numRuns); % LTRL2_LEC_V1
        
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
        
        % fmincon
%        s = s+1;
%         for ir = 1:numRuns
%            tStart = tic;
%            options_fmin_int.OutputFcn       = @(x,oV,state)(outputFcn_fmincon(x,oV,state,tStart));
%                                             
%                                         
%            [x_fmin,f_fmin,ex_fmin,out_fmin] = fmincon(@cutest_fun,x0,[],[],A,b,...,
%                                             [],[],[],options_fmin_int);
%                                         
%            time_alg         = toc(tStart); 
%            tcpu(s,ir,p)     = time_alg;
%            ex(p,s)          = ex_fmin;
%            numf(p,s)        = out_fmin.funcCount;
%            numg(p,s)        = out_fmin.funcCount;
%            numit(p,s)       = out_fmin.iterations;
%         end

        % Ipopt
        s = s+1;
        % exit : 0,1,2 acceptable
        for ir = 1:numRuns
          tStarti = tic;
          %funcs_ipopt.iterfunc       = @(x,oV,state)(...
          %                                  outputFcn_ipopt(x,oV,state,tStarti));
                                        
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
        
        p=p+1;            
    end
    
end

leg={   'TR1',...
        'TR1x1',...
        'TR1x2',...
        'TR2',...
        'TR2x1',...
        'TR2x2',...
        'IPOPT'};

% leg={   'TR1',...
%         'TR2',...       
%         'FMIN-LDL',...
%         'IPOPT'};

% leg={'TR-$(\mathbf{P},\infty)$',...
%        'TR-$\ell_2$',...
%        'fmincon-ldl'};
                
% types.colors    = ['b' 'r' 'm' 'g']; %'k' 'y'
% types.lines     = {'-', '-.', '-','-'}; %'-',   '-'
% types.markers   = ['o' 'o' 'o' 'o']; %'s' 's'
% 
% %indAlg          = [1 2 3 4 5]; %
% 
% %indAlg          = [1 2 3 4];
% indAlg          = [1 2 3];

types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

%indAlg          = [1 2 3 4 5]; %

indAlg          = [1 2 3 4 5 6 7]; % 5 6


perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_IV_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_IV_EXT.eps'));

% Comparison between proposed solvers only (TR1)
indAlg          = [1 2 3];

perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_IV_SEL_TR1_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_IV_SEL_TR1_EXT.eps'));

% Comparison between proposed solvers only (TR2)
indAlg          = [4 5 6];

perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_IV_SEL_TR2_EXT.eps'));

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
box on; grid on;
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_IV_SEL_TR2_EXT.eps'));

save(fullfile(datapath,'EXPERIMENT_IV_EXT'),'ex','numit','t_aver','numf','numg','params','tract', 'nms');

close ALL;

delete(     '*.d',...
            '*.o',...
            '*.dylib',...
            '*.f',...            
            'mcutest.*');
            
% Restore warning settings
warning(wtest);

tEX = toc(tEX); % Time to run experiment