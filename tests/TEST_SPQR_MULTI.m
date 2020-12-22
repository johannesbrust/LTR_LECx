%--------------------------- TEST_SPQR_MULTI -----------------------------%
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
optsSPQR.Q = 'Householder';

skipIdx = [65,66,70,71,0];
skipc = 1;
%fprintf('P# \t m      \t n      \t Time  \n');
fprintf('P# \t m      \t n      \t Time(F) \t Time(S)  \n');
for i = 1:nprob
    
    if i~=skipIdx(skipc)
    
        Prob= ssget(probIdx(i));
        A   = Prob.A;

        [m,n] = size(A);
        y = ones(n,1);

        tFAC = tic;
        [Q,~,~,info]=spqr(A',optsSPQR);
        tFAC = toc(tFAC);

        rkA         = info.rank_A_estimate;

        % Computing the projection via the Sparse QR
        t1      = tic;
        ztmp    = spqr_qmult(Q,y,0);
        ztmp1   = [zeros(rkA,1);ztmp(rkA+1:end)];
        z1      = spqr_qmult(Q,ztmp1,1);
        t1      = toc(t1);

        if m < 1000
            z2 = y - A'*pinv(full(A*A'))*(A*y);
            err = norm(z1-z2);
        end

        times(i,1) = tFAC;
        times(i,2) = t1;

        %fprintf('%i \t %i       \t %i       \t %.2e \n',i,m,n,t1);
        fprintf('%i \t %i       \t %i       \t %.2e \t %.2e \n',i,m,n,tFAC,t1);
    else
        skipc = skipc + 1;
    end
    
end
