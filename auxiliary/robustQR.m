function [ RA, maskA, rankA ] = robustQR( A, ranktol )
%ROBUSTQR Computes a robust LDLT factorization of A*A'
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
    [RA]    = qr(A);
else
    [RA]    = qr(A');
end

dRA     = diag(RA);

% Compute linearly independent columns.
maskA            = find(abs(dRA)>ranktol); % ^2
rankA            = length(maskA);  % rank of A

end

