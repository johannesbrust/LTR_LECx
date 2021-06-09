%%---------------------- PROBS_DENSECOLS --------------------------------%%
% From the article 'Large-Scale Quasi-Newton Trust-Region Methods
% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
% C.G. Petra
%
% This script checks for SuiteSparse problems with dense columns
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
%     06/08/21, J.B., Checking of dense columns

clc;
clear;

timeEX = tic;

addpath(genpath('../../main'));
addpath(genpath('../../auxiliary'));
%addpath ../solvers/rsqp
%addpath ../netlib/readmps
addpath(genpath('../../solvers/SPQR/MATLAB'));
addpath(genpath('../../ssget'));

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','..','/data/');
figpath     = fullfile(currentpath,'..','..','/figs/');
%probpath    = fullfile(currentpath,'..','..','/netlib/');

rng(090317);

fprintf('---------------- PROBS DENSECOLS ------------------------\n');

%%------------------ SuiteSparse Matrix Collection Problems --------------%
% Test matrices
%data = load('A_ARWHEAD');
data = load('EVEN_N_PROBS');
probIdx = data.probIdx;

%nprob = length(probIdx);
nprob = 50; % 30

msRks = data.msRks;

numn = nprob;%length(problems);
rks = zeros(numn,1);
nms = zeros(numn,2);
denseC = zeros(numn,1);

% Density tolerance
tolDense = 0.1;

optsSPQR.Q = 'Householder';

for i = 1:nprob
    % Loading problem and some data
    Prob= ssget(probIdx(i));
    A = Prob.A;
    
    [m,n] = size(A);
    
    % Computing a factorization for rank determination
    [Q,~,~,infoSP]=spqr(A',optsSPQR);
    rkA = infoSP.rank_A_estimate;
    
    rks(i) = rkA;
    nms(i,1) = m;
    nms(i,2) = n;
    
    on = true(n,1);
    om = ones(m,1);
    countD = (A'>0)*om;
    [vals,sIdx] = sort(countD,'descend');
    idxSi1 = sIdx(vals/m>tolDense); % Thresholding
    
    denseC(i) = length(idxSi1);
    
end

save('DENSE_COLS','rks','nms','denseC');

