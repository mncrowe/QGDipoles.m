function [K,a] = EVP(A,B,c,d,method,K_0)
% Solves the inhomogeneous EVP for n eignevalues:
%   (A - K*B])a = c_0 + K*c_1 s.t. d^T*a = 0
%
% A: matrix (N x N)
% B: matrix (N x N)
% c: vectors (N x 1), entered as array c = [c_0 c_1]
% d: vector (N x 1)
% method: 1 - invert B, 2 - invert A, 0 - invert matrix with larger rcond
% K_0: finds only the eigenvalue closest to K_0

% Invert either A or B, depending on which is better conditioned

if nargin < 5; method = 0; end

if (rcond(B) > rcond(A) || method == 1) && method ~= 2
    M1 = (d'*(B\c(:,1)))*A - (c(:,1)*d')*(B\A);
    M2 = -(d'*(B\c(:,1)))*B + (d'*(B\c(:,2)))*A - (c(:,2)*d')*(B\A);
    M3 = -(d'*(B\c(:,2)))*B;
    
else
    M1 = (d'*(A\c(:,1)))*A;
    M2 = (d'*(A\c(:,2)))*A - (d'*(A\c(:,1)))*B + (c(:,1)*d')*(A\B);
    M3 = -(d'*(A\c(:,2)))*B + (c(:,2)*d')*(A\B);
end

K = polyeig(M1,M2,M3);
K = K(~isinf(K));

if nargin > 5; [~,i] = min((K-K_0).^2); K = K(i); end

if nargout == 2
    a = zeros(sqrt(numel(A)),length(K));
    for i = 1:length(K)
        a(:,i) = (A-K(i)*B)\(c(:,1)+K(i)*c(:,2));
    end
end

end