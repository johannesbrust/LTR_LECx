function [ RA, maskA, rankA ] = robustLDLT( A, ranktol )
%ROBUSTLDLT Computes a robust LDLT factorization of A*A'
%
% Here A is (mxn) with m < n. 
%
% INPUTS:
% A : Rectancular matrix
% ranktol : Tolerance for rank determination
%
% OUTPUS:
% RA : Upper triangular factor, i.e., A*A' = RA'*RA
% maskA : Index of linearly independent rows of A
% rankA : Rank of A
%
%--------------------------------------------------------------------------
% 04/24/20, J.B., Initial implementation

[m,n] = size(A);

if m > n
    A = A';
    m = n;
end

sA      = sqrt(sum(A.^2,2));
idxnzr  = find(sA > ranktol);
%nnzr    = length(idxnzr);
%dsA     = diag(sA);
dsA     = spdiags(sA,0,m,m);
% Upper triangular matrix, and rank determination of A.
[Ldla, Ddla, pdla]  = ldl(dsA(idxnzr,(idxnzr))\...
                        ((A((idxnzr),:)*A((idxnzr),:)')/...
                        dsA((idxnzr),(idxnzr))),'vector');
                    
dDdla               = diag(Ddla);

% Compute linearly independent columns.
maskldla            = find(abs(dDdla)>ranktol); % ^2
rankA               = length(maskldla);  % rank of A

%pdlai(pdla(maskldla)) = 1:rankA;

maskA               = idxnzr(pdla(maskldla));  % index of safely linearly independent rows of A
%maskAi              = idxnzr(pdlai);

% Upper triangular factor A*A' = RA'*RA
RA(1:rankA,1:rankA)     = spdiags(sqrt(dDdla(maskldla)),0,rankA,rankA)*...
    (Ldla(maskldla,maskldla)'*dsA(maskA,maskA));


end

