%--------------------------- TEST_LSQR_MULTI -----------------------------%
%
% Testing the SuiteSparseQR methods from:
%
% https://github.com/DrTimothyAldenDavis/SuiteSparse/releases
%
% in order to compute projections z = P*y, where P = I - A'inv(A*A')A
%
% on a set of test problems from the SuiteSparseMatrix Collection. 
% This script references the "mex" files from SPQR/MATLAB
%-------------------------------------------------------------------------%
% 10/13/20, J.B., Initial test
% 10/15/20, J.B., Extension to use multiple A matrices
clc;
clear;

addpath(genpath('../solvers/SPQR/MATLAB'));
addpath('../ssget');

% Test matrices
%data = load('A_ARWHEAD');
data = load('EVEN_N_PROBS');
probIdx = data.probIdx;
%nprob = length(probIdx);
nprob = 60;

times = zeros(nprob,2);

% Options for SPQR
tolLS   = 1e-15;

optsLS.econ = 0;
optsLS.permutation = 'vector';
optsLS.Q = 'Householder';

fprintf('P# \t m      \t n      \t Time(F) \t Time(S)  \n');
for i = 1:nprob
    
    Prob= ssget(probIdx(i));
    A   = Prob.A;
    
    [m,n] = size(A);
    maxitLS = min(n,m);
    y = ones(n,1);
    
    
    tFAC = tic;
    [QLS,RLS,pLS,infoLS] = spqr(A',optsLS);
    tFAC = toc(tFAC);
    
    rkLS = infoLS.rank_A_estimate;
    linIND = pLS(1:rkLS);

    % Call LSQR
    tLS     = tic;
    [X,FLAG,RELRES,ITER] = lsqr(A(linIND,:)',y,tolLS,maxitLS,RLS(1:rkLS,1:rkLS)); %[]    
    zLS     = y - A(linIND,:)'*X;    
    tLS     = toc(tLS);
    
    times(i,1) = tFAC;
    times(i,2) = tLS;
    
    fprintf('%i \t %i       \t %i       \t %.2e \t %.2e \n',i,m,n,tFAC,tLS);
    
end
