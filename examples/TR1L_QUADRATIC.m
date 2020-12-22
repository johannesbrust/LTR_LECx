%-------------------------- TR1L_QUADRATIC -------------------------------%
%
% This example uses Algorithm 1 (with LSQR) {LTRL2_LEC_V2}
% on a quadratic objective:
%
% minimize c'*x + 0.5 x'*Q*x,   s.t. A*x = b.
%
% Here size(A) = [mm,n].
%
% Example for the algorithms from the report (or similar title): 
% 'Large-Scale Optimization with Linear Equality Constraints', 
% J.J. Brust, R.F. Marcia, C.G. Petra, M.A. Saunders, 2020
%
%-------------------------------------------------------------------------%
% 12/22/20, J.B.

clc;
clear;

addpath(genpath('../main'));
addpath(genpath('../auxiliary'));
addpath(genpath('../solvers/SPQR/MATLAB'));

n   = 1000;
mm  = 60; % 600

Q1  = randn(n,n);
Q   = Q1'*Q1;
A   = randn(mm,n);
b0  = randn(n,1);
b   = A*b0;
c   = randn(n,1);
x0  = A'*((A*A')\b);

fun         = @(x)( object_quad_arg(x,c,Q));
const       = @(x)( const_quad_arg(x,A,b));

options.storedat        = 0;
options.btol            = 1e-10;
options.gtol            = 1e-5;

options.trradb          = norm(x0);
    
options_l2              = options;
options_l2.maxitroot    = 10;
options_l2.epsroot      = 1e-5;
options_l2.dflag        = 1;

% Different possibilities of solving with A are:
%   whichAsolve - flag for which method to use in order to solve with A*A'
%                   1: LDLT factorization
%                   2: PCG iterative method
%                   3: QR factorization with Householder matrices
%                           (Uses SPQR from SuiteSparseQR)
%                   4: LSQR (with preconditioning)

% TR1L (using preconditioned LSQR)
options_l2.whichAsolve = 4;
[x_l2,f_l2,out_l2] = LTRL2_LEC_V2(fun,const,x0,options_l2);
