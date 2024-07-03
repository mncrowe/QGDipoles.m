% Solve for Beta-plane baroclinic topographic modons in a two layer model

clear; close all
addpath("functions","plotting")

% parameters:

U = 1;
a = 1;
R = [1 1];
beta = [0 1];
N = 2;
M = 12;

% derived quantities:

mu = beta*a^2/U;
lambda = a./R;

% create terms in eigenvalue problem:

A = zeros(N*M);
B = zeros(N*M);
Dn = @(n) kron(eye(M),diag(1:N == n));
c = [ones(N,1); zeros((M-1)*N,1)]/4;
dn = @(n) kron((-1).^(0:M-1)',(1:N == n)');

for j = 0:M-1
    for k = 0:M-1
        A(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) A_func(x,lambda,mu),j,k);
        B(j*N+(1:N),k*N+(1:N)) = JJ_int(@(x) B_func(x,lambda,mu),j,k);
    end
end

% solve eignevalue problem, (A-K_1^2*B_1-K_2^2*B_2)*a = (mu_1+K_1^2)*c_1 + 
% (mu_2+K_2^2)*c_2 s.t. d_1'*a = 0, d_2'*a = 0:

B1 = Dn(1)*B; B2 = Dn(2)*B;
Bcell = cell(2,1); Bcell{1} = B1; Bcell{2} = B2;
c_1 = Dn(1)*c; c_2 = Dn(2)*c;

[K,a_j] = EVP_optim(A,Bcell,[mu(1)*c_1+mu(2)*c_2 c_1 c_2],[dn(1) dn(2)],[4 5].^2,zeros(2*M,1));
K = sqrt(K);
a_j = reshape(a_j,[N M])';

% create plotting domain:

Nx = 512; Ny = 512;
Lx = 10; Ly = 10;
[x,y,k,l,grid] = create_domain([Nx Ny],[Lx Ly]);

% find psi_n and q_n using DeltaN(beta)*psi_n = -U/a*sin(theta)*sum[a_j R_j(r/a)] and q_n = DeltaN(0)*psi_n:

[psi_n,q_n] = CalcPsi(a_j,U,a,R,beta,grid);

% plot psi and q in each layer:

plot_window = [-2 2];
for n = 1:N
    Plot_2D(psi_n(:,:,n),x,y,plot_window,a)
    Plot_2D(q_n(:,:,n),x,y,plot_window,a)
end
