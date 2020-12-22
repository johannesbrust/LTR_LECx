function [ b ] = const_quad_arg_ipopt(x,A,b)
% CONST_QUAD_ARG_IPOPT constraint for a quadratic programming problem
% for use with Ipopt

    b  = A*x - b;

end