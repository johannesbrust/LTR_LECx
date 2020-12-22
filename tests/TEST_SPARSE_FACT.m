%--------------------- TEST_SPARSE_FACT ----------------------------------%
%
% Test for factorizing potentially large sparse rectangular systems.
%
% This test uses Matlab's documentation of 'Sparse Matrix Operations'.
% Including rank deficient cases.
%
% This is part of the extensions for the manuscript:
% 'Efficient Large-Scale Quasi-Newton Trust-Region Methods With 
% Linear Equality Constraints', J.J. Brust, R.F. Marcia, C.G. Petra, 2020
%-------------------------------------------------------------------------%
% 04/24/20, J.B.,
% 04/27/20, J.B., Modification of solution computation
% Representative times for n = 1e5, m = 2e3
%     Fact:                     Solve (Rankdef):
%     3.8590 (LDLT)             0.2504 (LDLT linsolve)
%     3.1387 (CHOL)             0.2718 (LDLT '\')
%    12.5704 (QR)               0.1651 (QR linsolve)
%     4.3702 (LDLT-Rankdef)     0.1778 (QR '\')
%    12.8745 (QR-Rankdef)       



addpath('../auxiliary')

n   = 1e5; % 1e6
m   = 2*1e3; % 1e5, 1*
rdef= 10;
mr  = m+rdef+1;
dens= 0.05; % 0.1
ntest = 6;
times= zeros(ntest,1);
ranktol = 1e-10;

A   = sprandn(m,n,dens);

% Rank deficient matrix
ARD = [A;...
        sprandn(rdef,m,dens)*A;...
        zeros(1,n)];

% Three tests: LDLT, chol, Q-less QR

solidx = 0;

%% LDLT factorization
solidx = solidx + 1;
ts = tic;
[RL,idxL,rankAL] = robustLDLT(A,ranktol);
times(solidx) = toc(ts); 

%% Cholesky factorization
solidx = solidx + 1;
ts = tic;
[RC,idxC,rankAC] = robustChol(A,ranktol);
times(solidx) = toc(ts); 

%% QR
solidx = solidx + 1;
ts = tic;
[RQ,idxQ,rankAQ] = robustQR(A,ranktol);
times(solidx) = toc(ts);

%% Rank deficient computations
%% LDLT factorization
solidx = solidx + 1;
ts = tic;
[RLr,idxLr,rankALr] = robustLDLT(ARD,ranktol);
times(solidx) = toc(ts); 

%% Cholesky factorization (Requires p.d. input)
% solidx = solidx + 1;
% ts = tic;
% [RCr,idxCr,rankACr] = robustChol(ARD,ranktol);
% times(solidx) = toc(ts); 

%% QR
solidx = solidx + 1;
ts = tic;
[RQr,idxQr,rankAQr] = robustQR(ARD,ranktol);
times(solidx) = toc(ts);

%% Solution times w. triangular matrices
opts1.UT        = true;
opts1.RECT      = false;
opts1.TRANSA    = false;
opts2           = opts1;
opts2.TRANSA    = true; % Transposed system
b               = randn(n,1);
solidx          = 0;
times2          = zeros(4,1);

%% ! Linsolve not supported for sparse inputs
% Linsolve LDLT
% Errors in solving rank deficient systems
solidx = solidx + 1;
ts = tic;

xlinL = ARD(idxLr,:)'*(linsolve(full(RLr(1:rankALr,1:rankALr)),...
     linsolve(full(RLr(1:rankALr,1:rankALr)),(ARD(idxLr,:)*b(:,1)),opts2),opts1));
times2(solidx) = toc(ts);

% Backslash LDLT
solidx = solidx + 1;
ts = tic;
%xBL = RL(1:rankAL,1:rankAL)\((RL(1:rankAL,1:rankAL)')\b(idxL));
xBL = ARD(idxLr,:)'*(RLr(1:rankALr,1:rankALr)\((RLr(1:rankALr,1:rankALr)')\(ARD(idxLr,:)*b(:,1))));
times2(solidx) = toc(ts);

% Linsolve QR
solidx = solidx + 1;
ts = tic;
xlinQ = ARD(idxQr,:)'*(linsolve(full(RQr(idxQr,idxQr)),...
     linsolve(full(RQr(idxQr,idxQr)),(ARD(idxQr,:)*b(:,1)),opts2),opts1));
times2(solidx) = toc(ts);

% Backslash QR
solidx = solidx + 1;
ts = tic;
%xBQ = RQ(idxQ,idxQ)\((RQ(idxQ,idxQ)')\b(idxQ));
xBQ = ARD(idxQr,:)'*(RQr(idxQr,idxQr)\((RQ(idxQr,idxQr)')\(ARD(idxQr,:)*b(:,1))));
times2(solidx) = toc(ts);

rIdxL(idxL) = 1:rankAL; 

err1 = norm(xlinL-xlinQ);
%err2 = norm(xBL(rIdxL)-xBQ);
err2 = norm(xBL-xBQ);
                                    













