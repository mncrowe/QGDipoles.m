function [x,y,k,l,grid] = create_domain(N,L)
% Creates Fourier grid
%
% N: [Nx Ny], number of gridpoints
% L: [Lx Ly], domain size

Nx = N(1); Ny = N(2);
Lx = L(1); Ly = L(2);

dx = Lx/Nx;
x = (-Lx/2:dx:(Lx/2-dx))';
k = (pi/dx*(-1:2/Nx:(1-2/Nx)))';

dy = Ly/Ny;
y = -Ly/2:dy:(Ly/2-dy);
l = pi/dy*(-1:2/Ny:(1-2/Ny));

if nargout > 4
    grid.x = x;
    grid.y = y;
    grid.k = k;
    grid.l = l;
end

end