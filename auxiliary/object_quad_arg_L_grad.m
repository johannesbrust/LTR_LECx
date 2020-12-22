function [ g ] = object_quad_arg_L_grad( phi, x, c, Q )
%OBJECT_QUAD_ARG_GRAD computes the function value and gradient
% of a quadratic objective function with symmetrix Hessian Q. For large
% problems.
%
%
%    The quadratic objective:
%
%    f = x'*c + 0.5* x'*Q*x,     g = c + Q*x 
%
% -----------------------------------------------------------------------
% 07/02/18, J.B., large version: Hessian = phi. I + Q Q'.
% 01/31/19, J.B., Gradient only

g = c + phi.*x + Q*(Q'*x);

end

