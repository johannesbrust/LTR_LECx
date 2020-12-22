function [ y ] = applyAA( x, A )
%APPLYAA computes A*A'*x
% Function to compute applies with rectangular A (mxn) with
% m < n. This is to be used with iterative solvers.
%
% INPUTS:
% x : Vector to be multiplied (mx1)
% A : Rectangular matrix (mxn)
%
% OUTPUTS:
% y : A*A'*x
%
%--------------------------------------------------------------------------
% 04/27/20, J.B., initial version

y = A*(A'*x);

end

