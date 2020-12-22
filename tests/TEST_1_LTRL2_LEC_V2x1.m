%------------------- TEST_1_LTRL2_LEC_V2x1 -------------------------------%
%
% Test 1 of the LTR2L2_LEC_V2 implementation with extension to use 
% SPQR to compute projections.
%
% This test is on small random possibly ('ill conditioned') quadratic 
% objective functions.
%
% This is part of the extensions for the manuscript:
% 'Efficient Large-Scale Quasi-Newton Trust-Region Methods With 
% Linear Equality Constraints', J.J. Brust, R.F. Marcia, C.G. Petra,
% M.A Saunders 2020
%-------------------------------------------------------------------------%
% 10/15/20, J.B.
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
options.dflag           = 0;
options.gtol            = 1e-5;

options.trradb          = norm(x0);
    
options_l2              = options;
options_l2.maxitroot    = 10;
options_l2.epsroot      = 1e-5;
options_l2.dflag        = 1;

% Version 2 (LSQR solver and projection based initialization)
options_l2.whichAsolve = 3;
options_l2.whichInit = 2;
[x_l2,f_l2,out_l2] = LTRL2_LEC_V2(fun,const,x0,options_l2);

options_l2.whichAsolve = 3;
options_l2.whichInit = 1;
[x_l2I1,f_l2I1,out_l2I1] = LTRL2_LEC_V2(fun,const,x0,options_l2);

%options_l2.whichAsolve = 1;
% Version 1
[x_l2v1,f_l2v1,out_l2v1] = LTRL2_LEC_V1(fun,const,x0,options_l2);