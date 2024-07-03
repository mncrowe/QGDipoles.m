function [psi_n, q_n] = CalcPsi(a_n,U,a,R,beta,grid)
% Calculate psi_n and q_n using D_N(beta)*psi = F, q = D_N(0)*psi.
%
% a_n: coefficients as M x N array
% (U,a): vortex parameters, scalars
% R: Rossby radii, vector of length N
% beta: beta, vector of length N
% grid: structure with fields x,y,k,l
%
% note: 'grid' is created using 5th output of 'create_domain.m'.

[M,N] = size(a_n);

% create RHS function F as a sum of Zernike radial functions:

r = sqrt(grid.x.^2+grid.y.^2);
theta = atan2(grid.y,grid.x);

F = zeros(length(grid.x),length(grid.y),N);
psi_n = zeros(length(grid.x),length(grid.y),N);
q_n = zeros(length(grid.x),length(grid.y),N);

for j = 1:M
    for i = 1:N
        F(:,:,i) = F(:,:,i) + a_n(j,i)*R_n(r/a,j-1);
    end
end
F = -U/a*F.*sin(theta);
F_FT = FT_2D(F,'forward');
F_FT(length(grid.x)/2+1,length(grid.y)/2+1,:) = 0;

% define operators D_N(beta)^-1 and D_N(0) in Fourier space:

D_inv = pageinv(DeltaN(grid.k,grid.l,R,beta/U));
D_q = DeltaN(grid.k,grid.l,R,0);

% calculate psi_n and q_n in Fourier space:

for n = 1:N
    for j = 1:N
        psi_n(:,:,n) = psi_n(:,:,n) + squeeze(D_inv(n,j,:,:)).*F_FT(:,:,j);
    end
end

for n = 1:N
    for j = 1:N
        q_n(:,:,n) = q_n(:,:,n) + squeeze(D_q(n,j,:,:)).*psi_n(:,:,j);
    end
end

% return psi_n and q_n to real space:

psi_n = real(FT_2D(psi_n,'inverse'));
q_n = real(FT_2D(q_n,'inverse'));

end