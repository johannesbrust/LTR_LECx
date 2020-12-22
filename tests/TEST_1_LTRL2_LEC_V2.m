%------------------ TEST_1_LTRL2_LEC_V2 ----------------------------------%
%
% Test 1 of the LTR2L2_LEC_V2 implementation.
%
% This test is on small random ('ill conditioned') quadratic objective
% functions.
%
% This is part of the extensions for the manuscript:
% 'Efficient Large-Scale Quasi-Newton Trust-Region Methods With 
% Linear Equality Constraints', J.J. Brust, R.F. Marcia, C.G. Petra, 2020
%-------------------------------------------------------------------------%
% 04/20/20, J.B.
clc;
clear;

addpath ../main
addpath ../auxiliary

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

% Version 2
[x_l2,f_l2,out_l2] = LTRL2_LEC_V2(fun,const,x0,options_l2);

% Version 1
[x_l2v1,f_l2v1,out_l2v1] = LTRL2_LEC_V1(fun,const,x0,options_l2);