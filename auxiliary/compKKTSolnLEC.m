function [ s ] = compKKTSolnLEC(delta,sigma,numsy,maskV,maskA,...
    nV,g,L,E,VV,V,A)
%COMPKKTSOLNLEC is is function to compute a KKT solution for
% testing purposes. Sigma enables the solution of shifted systems:
%
% [s; rho] = inv([ B + sigma.I, A'; A ])[-g; 0]
%
% INPUTS:
%  delta, sigma,..., and remaining inputs for computing B and 
%   bordered parts in the KKT system.
%
% OUTPUTS:
%   s:= Solution
%-------------------------------------------------------------------------%
% 04/21/20, J.B.

% Initializations
numsy2  = 2*numsy; 
n       = length(g);
mm      = length(maskA);
In      = eye(n);
zmm     = zeros(mm,mm);
M       = zeros(numsy2,numsy2);

% Explicit form of B
DV = diag(nV(maskV));

deltai = 1/delta;

M(1:numsy,1:numsy) = deltai.*VV(1:numsy,1:numsy);
M(1:numsy,(numsy+1):numsy2) = deltai.*L(1:numsy,1:numsy);
M((numsy+1):numsy2,1:numsy) = M(1:numsy,(numsy+1):numsy2)';
M((numsy+1):numsy2,(numsy+1):numsy2) = diag(-E(1:numsy));

M(1:numsy2,1:numsy2) = -M(1:numsy2,1:numsy2);

B = (delta+sigma).*In + ...
    V(:,maskV)*((DV(1:numsy2,1:numsy2)*M(1:numsy2,1:numsy2)*DV(1:numsy2,1:numsy2))\V(:,maskV)');

K = [B,A(maskA,:)';A(maskA,:), zmm(1:mm,1:mm)];
r = [-g;zmm(:,1)];

srho = K\r;

s    = srho(1:n);




