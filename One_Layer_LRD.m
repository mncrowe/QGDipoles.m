% Solve for the Larichev and Reznik Dipole (LRD) in a one layer QG model

clear; close all
addpath("functions","plotting")

% parameters:

U = 1;
a = 1;
R = 1;
beta = 1;
M = 12;

% derived quantities:

mu = beta*a^2/U;
lambda = a/R;

% create terms in eigenvalue problem:

A = zeros(M);
B = zeros(M);
c = [1/4; zeros(M-1,1)];
d = (-1).^(0:M-1)';

for j = 0:M-1
    for k = 0:M-1
        A(j+1,k+1) = JJ_int(@(x) A_func(x,lambda,mu),j,k);
        B(j+1,k+1) = JJ_int(@(x) B_func(x,lambda,mu),j,k);
    end
end

% solve eignevalue problem, (A-K^2*B)*a = (mu+K^2)*c s.t. d'*a = 0:

[K2,a2] = EVP(A,B,[mu*c c],d,1);

% extract the relevant solution, this is usually eigenvalue M-1:

Mode = M-1;
K = sqrt(K2(Mode)); a_j = a2(:,Mode);

% create plotting domain:

Nx = 512; Ny = 512;
Lx = 10; Ly = 10;
[x,y,k,l,grid] = create_domain([Nx Ny],[Lx Ly]);

% find psi_n and q_n using DeltaN(beta)*psi_n = -U/a*sin(theta)*sum[a_j R_j(r/a)] and q_n = DeltaN(0)*psi_n:

[psi_n,q_n] = CalcPsi(a_j,U,a,R,beta,grid);

% plot psi and q in each layer:

plot_window = [-2 2];
Plot_2D(psi_n,x,y,plot_window,a)
Plot_2D(q_n,x,y,plot_window,a)
