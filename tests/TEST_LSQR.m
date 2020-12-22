%--------------------------- TEST_LSQR -----------------------------------%
%
% Testing the LSQR method for integration into TR1, TR2
%
% This script references the "mex" files from SPQR/MATLAB
%-------------------------------------------------------------------------%
% 10/14/20, J.B., Initial test

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

% Computation using LSQR and preconditioning
tolLS   = 1e-15;
maxitLS = min(n,m);

optsLS.econ = 0;
optsLS.permutation = 'vector';

%RLS     = spqr(A',optsLS);
[QLS,RLS,pLS,infoLS] = spqr(A',optsLS);
rkLS = info.rank_A_estimate;

linIND = pLS(1:rkLS);

% Set nearly zero diagonal to 1
%idxLrk  = find(abs(diag(RLS))<1e-10);
% idxLrk  = find(abs(diag(RLS))>1e-10);
% idxD    = sub2ind(size(RLS),idxLrk,idxLrk);
%RLS(idxD) = 1;

% Call LSQR
tLS     = tic;
%[X,FLAG,RELRES,ITER] = lsqr(A(idxLrk,:)',y,tolLS,maxitLS,RLS(idxLrk,idxLrk));
[X,FLAG,RELRES,ITER] = lsqr(A(linIND,:)',y,tolLS,maxitLS,RLS(1:rkLS,1:rkLS)); %[]
%[X,FLAG,RELRES,ITER] = lsqr(A',y,tolLS,maxitLS,[]); %[]
%zLS     = y - A(idxLrk,:)'*X;
zLS     = y - A(linIND,:)'*X;
%zLS     = y - A'*X;
tLS     = toc(tLS);

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
err2    = norm(zLS-z2);