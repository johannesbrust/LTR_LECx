%--------------------- TEST_SPARSE_ITSOLVE -------------------------------%
%
% Test for solving with large rectangular linear systems. In addition
% to a LDLT factorization approach, the iterative solvers 'pcg', 'minres'
% and 'symmlq' are used.
%
% This test uses Matlab's documentation of 'Sparse Matrix Operations'.
% Including rank deficient cases.
%
% This is part of the extensions for the manuscript:
% 'Efficient Large-Scale Quasi-Newton Trust-Region Methods With 
% Linear Equality Constraints', J.J. Brust, R.F. Marcia, C.G. Petra, 2020
%-------------------------------------------------------------------------%
% 04/27/20, J.B.,
% Representative result: n = 1e6, m = 1e4
% Meth      Full rank:     Rank deficient:
% PCG           5.4129         13.3730
% MINR          5.4092         15.7644
% SYMM          5.7735         15.4586

addpath('../auxiliary')
clc;
clear;

tdatas = tic;
n   = 1e5; % 1e6, 1e5, 1e6
m   = 1e3; % 1e5, 3*1e3, 1e4
rdef= 10;
mr  = m+rdef+1;
dens= 0.02; % 0.1
ntest = 4;
times= zeros(ntest,2); % First column full-rank, second column not full rank
timesfac = zeros(2,1); % Factorization times
ranktol = 1e-10;

A   = sprandn(m,n,dens);

% Rank deficient matrix
ARD = [A;...
        sprandn(rdef,m,dens)*A;...
        zeros(1,n)];

% Right hand side and function handles
g       = randn(n,1);
g(n)    = 0;
aFunc   = @(x)(applyAA(x,A));
ardFunc = @(x)(applyAA(x,ARD));

Ag      = A(:,:)*g(:);
ARDg    = ARD(:,:)*g(:);
nAg     = norm(Ag);
nARDg   = norm(ARDg);

maxit   = 500;
tolA    = max(1e-9/nAg,1e-12);
tolARD  = max(1e-9/nARDg,1e-12);
xinitA  = zeros(m,1);
xinitARD= zeros(mr,1);

tdata = toc(tdatas);
% Four methods: LDLT, pcg, minres, symmlq

solidx = 0;

%% LDLT factorization and solve
if m*n < 5e8
    solidx = solidx + 1;
    tfacs   = tic;
    ts      = tic;
    % Factorization
    [RL,idxL,rankAL] = robustLDLT(A,ranktol);
    timesfac(1) = toc(tfacs);
    % Solve
    xBL = A(idxL,:)'*(RL(1:rankAL,1:rankAL)\((RL(1:rankAL,1:rankAL)')\(A(idxL,:)*g(:))));
    times(solidx,1) = toc(ts); 
end

%% PCG
solidx = solidx + 1;
ts = tic;
Ag = A(:,:)*g(:);
[xPCG,fPCG,resPCG,iterPCG] = pcg(aFunc,Ag,tolA,maxit,[],[],xinitA);
xPCGo(:,1) = A(:,:)'*xPCG(:,1);
times(solidx,1) = toc(ts); 

%% MINRES
solidx = solidx + 1;
ts = tic;
Ag = A(:,:)*g(:);
[xMINR,fMINR,resMINR,iterMINR] = minres(aFunc,Ag,tolA,maxit,[],[],xinitA);
xMINRo(:,1) = A(:,:)'*xMINR(:,1);
times(solidx,1) = toc(ts); 

%% SYMMLQ
solidx = solidx + 1;
ts = tic;
Ag = A(:,:)*g(:);
[xSYM,fSYM,resSYM,iterSYM] = symmlq(aFunc,Ag,tolA,maxit,[],[],xinitA);
xSYMo(:,1) = A(:,:)'*xSYM(:,1);
times(solidx,1) = toc(ts); 

%% Rank deficient computations
solidx = 0;
%% LDLT factorization
if m*n < 5e8
    solidx  = solidx + 1;
    tfacs   = tic;
    ts      = tic;
    [RLr,idxLr,rankALr] = robustLDLT(ARD,ranktol);
    timesfac(2) = toc(tfacs);
    % Solve
    xBLr = ARD(idxLr,:)'*(RLr(1:rankALr,1:rankALr)\((RLr(1:rankALr,1:rankALr)')\(ARD(idxLr,:)*g(:))));
    times(solidx,2) = toc(ts); 
end

%% PCG
solidx = solidx + 1;
ts = tic;
ARDg = ARD(:,:)*g(:);
[xPCGr,fPCGr,resPCGr,iterPCGr] = pcg(ardFunc,ARDg,tolARD,maxit,[],[],xinitARD);
xPCGor(:,1) = ARD(:,:)'*xPCGr(:,1);
times(solidx,2) = toc(ts); 

%% MINRES
solidx = solidx + 1;
ts = tic;
ARDg = ARD(:,:)*g(:);
[xMINRr,fMINRr,resMINRr,iterMINRr] = minres(ardFunc,ARDg,tolARD,maxit,[],[],xinitARD);
xMINRor(:,1) = ARD(:,:)'*xMINRr(:,1);
times(solidx,2) = toc(ts); 

%% SYMMLQ
solidx = solidx + 1;
ts = tic;
ARDg = ARD(:,:)*g(:);
[xSYMr,fSYMr,resSYMr,iterSYMr] = symmlq(ardFunc,ARDg,tolARD,maxit,[],[],xinitARD);
xSYMor(:,1) = ARD(:,:)'*xSYMr(:,1);
times(solidx,2) = toc(ts); 



% 
% %% Solution times w. triangular matrices
% opts1.UT        = true;
% opts1.RECT      = false;
% opts1.TRANSA    = false;
% opts2           = opts1;
% opts2.TRANSA    = true; % Transposed system
% b               = randn(n,1);
% solidx          = 0;
% times2          = zeros(4,1);
% 
% %% ! Linsolve not supported for sparse inputs
% % Linsolve LDLT
% solidx = solidx + 1;
% ts = tic;
% 
% xlinL = ARD(idxLr,:)'*(linsolve(full(RLr(1:rankALr,1:rankALr)),...
%      linsolve(full(RLr(1:rankALr,1:rankALr)),(ARD(idxLr,:)*b(:,1)),opts2),opts1));
% times2(solidx) = toc(ts);
% 

                                    













