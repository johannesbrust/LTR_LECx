%%----------------------- TEST_EXT_1a -------------------------------------
%% From the article 'Large-Scale Quasi-Newton Trust-Region Methods
%% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
%% C.G. Petra. EXTENSION: Large-scale linear equality constraints.
%
% This test implements a formula for the search direction that requires
% only one solve with AA' per iteration.
% 
% TEST_EXT_1a implements correction terms based on the orthogonal
% projection
%
% skp = (I-A'inv(A*A')A)*sk
%
%--------------------------------------------------------------------------
% 03/18/20, J.B., Initial implementation
% 03/19/20, J.B., Test including the L2 solver
%               , Correction term
% 03/25/20, J.B., Correction using an orthogonal projection
%               , KKT step

clc;
clear;

addpath ../main
addpath ../auxiliary

% Initializations
n = 100;
m = 50;
maxiter = 20;
epscorr = 1e-8;

nerr = zeros(maxiter+1,1);

% Problem data
Q1      = randn(n,n);
Q       = Q1'*Q1;
A       = randn(m,n);
xk      = 5.*ones(n,1);
x0      = xk;
b       = A*xk;
gk      = Q*xk;
Agk     = A*gk;
AAigk   = (A*A')\Agk;

% Algorithm storage
S       = zeros(n,maxiter);
Y       = zeros(n,maxiter);
E       = zeros(maxiter,maxiter);
YY      = zeros(maxiter,maxiter);
T       = zeros(maxiter,maxiter);
%AY      = zeros(m,maxiter);
AAiAY   = zeros(m,maxiter);
YAAAiAY = zeros(maxiter,maxiter);
M       = zeros(2*maxiter,2*maxiter);

L       = zeros(maxiter,maxiter);

sk      = -(gk - A'*AAigk);

xk1     = xk + sk;
gk1     = Q*xk1;
rk      = A*sk;
%rk      = A*xk1-b;

yk      = gk1-gk;
Agk1    = A*gk1;
AAigk1  = (A*A')\Agk1;

delta   = 1;

Imax    = eye(maxiter);
In      = eye(n);
zm      = zeros(m,m); 

nrk     = norm(rk);
nerr(1) = nrk;
correct = 0;

B       = In;
solK    = [B, A'; A, zm ]\[-gk;zm(:,1)];

fprintf('EXTENDED LTR METHODS \n');
fprintf('Iter \t norm(rk)  \t norm(gk) \t norm(sk) \t diff(sk) \n');
fprintf('%i \t %.3g \t %.3g \t %.3g \t %.3g \n',1,nrk,norm(gk1),norm(sk),norm(sk-solK(1:n)));

% Main iteration
for k = 1:maxiter
   
    % Updates
    S(:,k)  = sk;
    Y(:,k)  = yk;
    AAiAY(:,k) = AAigk1 - AAigk;
    YAAAiAY(1:k,1:k) = (Y(:,1:k)'*A')*AAiAY(:,1:k);
    
    SY          = S(:,1:k)'*Y(:,1:k);
    
    T(1:k,1:k)  = triu(SY,0);
    L(1:k,1:k)  = SY - T(1:k,1:k);
    E(1:k,1:k)  = diag(diag(T(1:k,1:k)));
    
    % Step computation
    gk      = gk1;
    xk      = xk1;
    Agk     = Agk1;
    AAigk   = AAigk1;
    
    SYgk    = [S(:,1:k)'*gk;Y(:,1:k)'*gk];
    
    p11 = -AAigk;
    
    M   = [((T(1:k,1:k)')\(E(1:k,1:k)+Y(:,1:k)'*Y(:,1:k)))/T(1:k,1:k), -(T(1:k,1:k)\Imax(1:k,1:k))';...
            -(T(1:k,1:k)\Imax(1:k,1:k)), zeros(k,k)];
        
    p12 = -AAiAY(:,1:k)*(M((k+1):2*k,1:2*k)*SYgk);
    
    p21 = - M(:,(k+1:2*k))*(AAiAY(:,1:k)'*Agk);
    
    M1 = [((T(1:k,1:k)')\(E(1:k,1:k)+Y(:,1:k)'*Y(:,1:k)-YAAAiAY(1:k,1:k)))/T(1:k,1:k), -(T(1:k,1:k)\Imax(1:k,1:k))';...
            -(T(1:k,1:k)\Imax(1:k,1:k)), zeros(k,k)];
    
    p22 = M1*SYgk;
        
    
    sk  = -A'*(p11 + p12) - [S(:,1:k), Y(:,1:k)]*(p21 + p22) - gk;
    
    % Corrections based on an orthogonal projection
    
    %cr1 = A'*AAigk;
    %cr3 = A'*(AAiAY(:,1:k)*(T(1:k,1:k)\(-SYgk(1:k))));
    
    %sk  = sk_ + cr1 + cr3;
    
    %sk   = sk_ -A'*((A*A')\(A*sk_));
    
    % KKT step
    %B       = In - [S(:,1:k), Y(:,1:k)]*([S(:,1:k)'*S(:,1:k), L(1:k,1:k);...
    %                L(1:k,1:k)', -E(1:k,1:k)]\([S(:,1:k), Y(:,1:k)]'));
                
    %solK    = [B, A'; A, zm ]\[-gk;zm(:,1)];            
    %sk = solK(1:n);
    
    alp = -(sk'*gk)/(sk'*Q*sk);
    
    xk1 = xk + alp.*sk;
    %xk1 = xk + solK(1:n);
    
    rk  = A*sk;
    %rk  = A*xk1-b;
    nrk = norm(rk);
    
    gk1 = Q*xk1;

    yk  = gk1 - gk;
    
    Agk1    = A*gk1;
    AAigk1  = (A*A')\Agk1;
    
    nerr(k+1) = nrk;
    %fprintf('%i \t %.3g \t %i \n',k+1,nrk,correct);    
    fprintf('%i \t %.3g \t %.3g \t %.3g \t %.3g \n',k+1,nrk,norm(gk1),norm(sk),norm(sk-solK(1:n)));

    correct = 0;
    
end

em = A*S(:,1:maxiter);

%% L2 trust-region method

options_const.storedat  = 0;
options_const.btol      = 1e-10;
options_const.dflag     = 0;
options_const.gtol      = 5e-5;
options_const.maxitroot = 10;
options_const.epsroot   = 1e-5;

c                       = zeros(n,1);

obj                     = @(x)( object_quad_arg(x,c,Q));
const                   = @(x)( const_quad_arg(x,A,b));

[xl2,fl2,outinfol2]     = LTRL2_LEC_V1(obj,const,x0,options_const);



































