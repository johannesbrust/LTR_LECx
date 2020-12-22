function [ g ] = object_quad_arg_grad( x, c, Q )
%OBJECT_QUAD_ARG_GRAD computes the function gradient
% of a quadratic objective function with symmetrix Hessian Q. For use with
% Ipopt

%{

    06/21/18 J.B.

    f = x'*c + 0.5* x'*Q*x,     g = c + Q*x 

%}

g = c + Q*x;


end

