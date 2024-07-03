% solve for a general N layer modon
%
% This script is currently set up for a 5 layer modon. This can changed by
% modifying R and beta. Length(R), length(beta) = N. When entering K_0, use
% values slightly larger than you expect for K, this seems to converge. The
% most likely issue here is the solution having a passive layer where an
% active one is requested. This generally means the initial K_0 is not
% close enough to the solution for K with an active layer. Try increasing
% the corresponding value of K_0.
%
% Note: Instead of using this general approach, it may be easier to find
% solutions with passive layers by explicitly removing those K_i and only
% solving for the active layer eigenvalues. However, the EVP will still
% have the same number of parameters to solve for as each extra K_i also
% introduces a condition which reduces the number of free parameters in a.

clear; close all
addpath("functions","plotting")

% physical parameters:

U = 1;                      % vortex speed
a = 1;                      % vortex radius
R = [1 1 2 1 2];            % Rossby radius in each layer
beta = [1 0 1 2 0];         % background vorticity gradients
AL = [1 1 0 1 1];           % layer type; 1-active, 0-passive

% numerical parameters:

M = 5;                     % number of coefficients in a_n
K_0 = [5 5 4 5 7];          % guess for K in each layer
show_diag = 0;              % 1 - show fsolve output for EVP

% derived quantities:

N = length(R);
mu = beta*a^2/U;
lambda = a./R;
K_0(~AL) = 1i*sqrt(mu(~AL));    % set -K_i^2 = mu_i in passive layer

% create terms in eigenvalue problem:

A = zeros(N*M);
B = zeros(N*M);
c = zeros(N*M,N+1);
d = zeros(N*M,N);
Bcell = cell(N,1);

for j = 0:M-1
    for k = 0:M-1
        A(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) A_func(x,lambda,mu),j,k);
        B(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) B_func(x,lambda,mu),j,k);
    end
end

for n = 1:N
    Bcell{n} = kron(eye(M),diag(1:N == n))*B;
    c(:,1) = c(:,1) + mu(n)*kron(eye(M),diag(1:N == n))*[ones(N,1); zeros((M-1)*N,1)]/4;
    c(:,n+1) = kron(eye(M),diag(1:N == n))*[ones(N,1); zeros((M-1)*N,1)]/4;
    d(:,n) = kron((-1).^(0:M-1)',(1:N == n)');
end

% solve eignevalue problem with initial guess K = K_0:

[K,a_n,err] = EVP_optim(A,Bcell,c,d,K_0.^2,zeros(N*M,1),show_diag);
K = sqrt(K);
a_n = reshape(a_n,[N M])';

% display solution details:

AL2 = abs(K.^2+mu)>1e-2;
disp(' ')
disp(['Solution found for (U,a) = (' num2str(U) ',' num2str(a) ') with N = ' num2str(N) ' layers and M = ' num2str(M) ' coefficients:'])
disp(' ')
for n = 1:N
    if AL2(n) == 1
        disp([' - Layer ' num2str(n) ' with (R,beta) = (' num2str(R(n)) ',' num2str(beta(n)) '): Active with K_i = ' num2str(K(n))])
    else
        disp([' - Layer ' num2str(n) ' with (R,beta) = (' num2str(R(n)) ',' num2str(beta(n)) '): Passive with mu_i = ' num2str(mu(n))])
    end
end
if ~isequal(AL2,AL); disp(' '); disp('Active/passive layers do not match requested solution, try a different K_0.'); end
disp(' ')

% create plotting domain:

Nx = 512; Ny = 512;
Lx = 20; Ly = 20;
[x,y,k,l,grid] = create_domain([Nx Ny],[Lx Ly]);

% find psi_n, q_n using DeltaN(beta)*psi_n = -U/a*sin(theta)*sum[a_j R_j(r/a)], q_n = DeltaN(0)*psi_n:

[psi_i,q_i] = CalcPsi(a_n,U,a,R,beta,grid);

% plot psi and q in each layer:

plot_window = [-2 2];
for n = 1:N
    Plot_2D(psi_i(:,:,n),x,y,plot_window,a)
    Plot_2D(q_i(:,:,n),x,y,plot_window,a)
end
