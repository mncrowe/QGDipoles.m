function [K,a,err] = EVP_optim(A,B,c,d,K_0,a_0,diag,tol)
% Solves the problem (A-sum[K_i*B_i])*a = c_0 + sum[K_i*c_i] such that
% d_k^T*a = 0 for i,k = {1,2,...,N}.
%
% A: M x M matrix
% B: cell array of N M x M matrices (B = cell(N,1); B{i} = B_i)
% c: vectors (M x 1), enter as array c = [c_0 c_1 c_2 .. c_N]
% d: vectors (M x 1), enter as array d = [d_1 d_2 .. d_N]
% K_0: initial guess for [K_1, K_2, ..., K_N]
% a_0: initial guess for a, vector of length M
% diag: shows iteration output if set to 1 (default: 0)
% tol: function tolerance for fsolve (default: 1e-10)
%
% This function works by defining:
% 
%       f(x) = (A-sum[K_i*B_i])*a - c_0 - sum[K_i*c_i],
%
% where x = [K_1, K_2, .. K_N, a'_1, .. a'_{M-N}] and the a'_i are
% coefficients of the M-N basis vectors which span the space orthogonal to
% [d_1, d_2, ..., d_N]. We now find the solution to f(x) = 0  using the
% Matlab fsolve function from the Optimisation tooolbox. The solution x
% gives the values of K_i and allows a to be reconstructed from the a'_i. 
% Note that the complexity of this problem does not increase as the number
% of eigenvalues, N, does since a new eigenvalue introduces a new condition
% which reduces the dimension of the space containing a. Resonable guesses
% for K_0 are generally sufficient, allowing a_0 to be left as 0.

% first find a basis for the space orthogonal to the d_j:

[e,V,I] = OrthogSpace(d);
[N,J] = size(e);

% create initial guess randomly or from supplied data if given:

if nargin == 5; a_0 = zeros(N,1); end
if nargin > 4
    x0 = V\a_0;
    x0 = [reshape(K_0,N-J,1); x0(I)];
else
    x0 = rand(N,1);
end
if nargin < 7; diag = 0; end
if nargin < 8; tol = 1e-10; end

% find solution x which solves f(x) = 0:

f = @(x) EVP_optim_f(x,A,B,c,e);
if diag == 1
    options = optimoptions('fsolve','Display','iter','TolFun',tol,'SpecifyObjectiveGradient',true);
else
    options = optimoptions('fsolve','Display','none','TolFun',tol,'SpecifyObjectiveGradient',true);
end
[x,err] = fsolve(f,x0,options);

% determine K and a from x:

K = reshape(x(1:N-J),1,N-J);
a = e*x(N-J+1:N);

end

% definition of the function f(x):

function [f_out,Jac] = EVP_optim_f(x,A,B,c,e)
       
    [N,J] = size(e);

    a = e*x(N-J+1:N);
    M = A;
    v = c(:,1);
    Jac = zeros(N,N-J);

    for i = 1:N-J
        M = M - x(i)*B{i};
        v = v + x(i)*c(:,i+1);
        Jac(:,i) = -B{i}*a-c(:,i+1);
    end

    f_out = (M*a-v);
    Jac = [Jac M*e];
    
end