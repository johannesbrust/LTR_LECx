%--------------------------- TEST_SPQR -----------------------------------%
%
% Testing the SuiteSparseQR methods from:
%
% https://github.com/DrTimothyAldenDavis/SuiteSparse/releases
%
% in order to compute projections z = P*y, where P = I - A'inv(A*A')A.
%
% This script references the "mex" files from SPQR/MATLAB
%-------------------------------------------------------------------------%
% 10/13/20, J.B., Initial test

addpath(genpath('/Users/johannesbrust/Dropbox/codes/SuiteSparse-5.8.1/SPQR/MATLAB'));

% Test matrix and vector
rng('default');
m       = 100;
n       = 2000;
dens    = 0.05;
A       = sprandn(m,n,dens);
y       = ones(n,1);

rkDef   = 1;
if rkDef == 1
    A(10,:) = zeros(1,n);
    A(99,:) = 2*A(100,:);
end

% Options for SPQR
optsSPQR.Q = 'Householder';

[Q,~,~,info]=spqr(A',optsSPQR);
rkA         = info.rank_A_estimate;

% Computing the projection via the Sparse QR
t1      = tic;
ztmp    = spqr_qmult(Q,y,0);
ztmp1   = [zeros(rkA,1);ztmp(rkA+1:end)];
z1      = spqr_qmult(Q,ztmp1,1);
t1      = toc(t1);

% Computing the projection directly
t2      = tic;
z2      = y - A'*(pinv(full(A*A'))*(A*y));
t2      = toc(t2);

err1    = norm(z1-z2);