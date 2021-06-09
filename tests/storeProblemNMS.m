%-------------------- storeProblemNMS ------------------------------------%
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
% This script stores the problem dimensions (m,n)
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
% 06/04/21, J.B., Storing problem data
%--------------------------------------------------------------------------
%NOTE(S): This Experiment requires CUTEst to be installed

clc;
clear;

addpath(genpath('../main'));
addpath(genpath('../auxiliary'));
addpath(genpath('../netlib/readmps'));
addpath(genpath('../solvers/SPQR/MATLAB'));

if ismac == 1
    addpath(genpath('../solvers/IPOPT_UP'));
else
    addpath(genpath('../solvers/ipopt'));
end

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','/data/');
figpath     = fullfile(currentpath,'..','/figs/');
probpath    = fullfile(currentpath,'..','/auxiliary/');
%probpath    = fullfile(currentpath);

CUTEst_init  % initialize CUTEr, see appropriate documentation
%fid     = fopen(fullfile(probpath,'probs_EX_V.txt'),'r');
fid     = fopen(fullfile(probpath,'cutest_list_EX_V.txt'),'r');

numProblems = 31;
nms             = zeros(numProblems,2);
p=1;

tline = fgets(fid);

while ischar(tline)
        
    if (~strcmp(tline(1),'%'))  && (ischar(tline))
        
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
        n               = size(x0,1);
        
        % CUTEST constraint
        [c,Ad] = cutest_cons(x0);
        A = sparse(Ad);
        mm = size(A,1);        
        
        nms(p,1) = n;
        nms(p,2) = mm;
        
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

save('problemNMS','nms');