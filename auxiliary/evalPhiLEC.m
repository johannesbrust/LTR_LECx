function [ phis, dphis ] = evalPhiLEC(sigs,numsy,m,delta,VV,T,E,Vg,V,gp,nV,...
    maskV,opts3,trrad)
%EVALPHILEC is a function to evaluate the function
% phi(sig) = 1/nsk - 1/Deltak in the TR subproblem for
% testing purposes.
%
% INPUTS:
%   sigs:=      Sigma values to evaluate phi
%
%   numsy,..., and remaining inputs for computing phi and its derivative.
%
% OUTPUTS:
%   phis:=      Phi values
%   dphis:=     Computed Phi' values
%-------------------------------------------------------------------------%
% 04/14/20, J.B., Initial Implementation

numsy2  = 2*numsy;

nsigs   = length(sigs);

phis    = zeros(nsigs,1);
dphis   = zeros(nsigs,1);

N       = zeros(2*m,2*m);

for i = 1:nsigs
    
    sigma = sigs(i);
    
    gam     = delta+sigma;
    del     = 1/gam;
    alp     = delta/gam;
    a1      = 1/alp-1;
    
    %N(1:numsy,1:numsy)
    %     N(1:numsy,(m+1):(m+numsy))  = a1.*VV(1:numsy,m+1:m+numsy)-T(1:numsy,1:numsy)./alp;
    %     N((m+1):(m+numsy),1:numsy)  = N(1:numsy,(m+1):(m+numsy));
    %     N(1:numsy,1:numsy)          = a1.*VV(1:numsy,1:numsy);
    %     N(m+1:m+numsy,m+1:m+numsy)  = -(diag(gam.*E(1:numsy))+VV(m+1:m+numsy,m+1:m+numsy));
    
    N(1:numsy,(numsy+1):numsy2)  = a1.*VV(1:numsy,m+1:m+numsy)-T(1:numsy,1:numsy)./alp;
    N((numsy+1):numsy2,1:numsy)  = N(1:numsy,(numsy+1):numsy2)';
    N(1:numsy,1:numsy)          = a1.*VV(1:numsy,1:numsy);
    N((numsy+1):numsy2,(numsy+1):numsy2)  = -(diag(gam.*E(1:numsy))+VV(m+1:m+numsy,m+1:m+numsy));
    
    p(1:numsy2,1)  = linsolve(N(1:numsy2,1:numsy2),Vg(1:numsy2),opts3);
    
    st           = -del.*(V(:,maskV)*(p(1:numsy2)./nV(maskV)) + gp);
    
    Vs(1:numsy2,1)            = (V(:,maskV)'*st)./nV(maskV);
    ps(1:numsy2,1)  = linsolve(N(1:numsy2,1:numsy2),Vs(1:numsy2,1),opts3);
    
    nst                 = norm(st);
    
    %     nst_                = del*sqrt(ngp^2+2*(p(1:numsy2,1)'*Vg(1:numsy2))+...
    %         p(1:numsy2,1)'*(VV([1:numsy,m+1:m+numsy],[1:numsy,m+1:m+numsy])*p(1:numsy2,1)));
    
    phi                 = 1/nst - 1/trrad;
    
    stpst               = -del*(Vs(1:numsy2)'*ps(1:numsy2,1) + nst^2);
    
    phip                =-stpst/nst^3;
    
    phis(i) = phi;
    dphis(i) = phip;
    
end


