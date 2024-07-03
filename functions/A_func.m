function A = A_func(xi,lambda,mu)
% Calculates the integrand for A_kj, K(xi)*(K(xi)+D(mu))^-1*xi^-1
%
% xi: coordinate vector
% lambda: column vector of lambda values
% mu: column vector of mu values

lambda = reshape(lambda,[],1);

N = length(mu);
if N == 1
    K = xi.^2+lambda^2;
    A = K./(K+mu)./xi;
else
    K = repmat(eye(N),1,1,length(xi)).*reshape(xi.^2,1,1,length(xi))+...
        diag([lambda(1)^2 2*reshape(lambda(2:end-1).^2,1,[]) lambda(end)^2])+...
        diag(-lambda(1:end-1).^2,1)+diag(-lambda(2:end).^2,-1);
    A = pagemrdivide(K,K+diag(mu))./reshape(xi,1,1,length(xi));
end

end