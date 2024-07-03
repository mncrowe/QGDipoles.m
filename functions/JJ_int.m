function I = JJ_int(F,j,k,tol)
% Performs the double Bessel function integral int_0^inf F(z) J_{2j+2}(z) J_{2k+2}(z) dz
%
% F: function
% j,k: Bessel function orders
% tol: error tolerance, default: 1e-6
%
% Domain is [0,inf] and is split at d. For x < d, the contribution is
% evalulated directly, for x > d we use a contour deformation approach to
% avoid strongly oscillatory behaviour and slow convergence. Method
% suggested by A. Gibbs.

d = 1000;               % domain splitting parameter
N = sqrt(numel(F(0)));  % dimension of F (i.e. F : R^{N x N} -> R)
if nargin < 4; tol = 1e-6; end
atol = 1e-4*tol;

if N == 1
    R = @(x) x;
else
    R = @(x) reshape(x,1,1,[]);
end

J = @(n,x) besselj(n,x);
Y = @(n,x) bessely(n,x);
H1 = @(n,x) besselh(n,1,x);
H2 = @(n,x) besselh(n,2,x);
HH12 = @(m,n,x) 2*(J(m,x).*J(n,x)+Y(m,x).*Y(n,x));
phi11 = @(m,n,x) H1(m,x).*H1(n,x)./exp(2*1i*x);
phi22 = @(m,n,x) H2(m,x).*H2(n,x)./exp(-2*1i*x);

I1 = integral(@(x) F(x).*R(J(2*j+2,x).*J(2*k+2,x)),0,d,'ArrayValued',true,'AbsTol',atol,'RelTol',tol);
I2 = integral(@(x) F(x).*R(HH12(2*j+2,2*k+2,x)/4),d,inf,'ArrayValued',true,'AbsTol',atol,'RelTol',tol);
I3 = 1i*exp(2*1i*d)/4*integral(@(x) F(d+1i*x).*R(phi11(2*j+2,2*k+2,d+1i*x).*exp(-2*x)),0,100,'ArrayValued',true,'AbsTol',atol,'RelTol',tol);
I4 = -1i*exp(-2*1i*d)/4*integral(@(x) F(d-1i*x).*R(phi22(2*j+2,2*k+2,d-1i*x).*exp(-2*x)),0,100,'ArrayValued',true,'AbsTol',atol,'RelTol',tol);

I = I1+I2+real(I3+I4);

end