function F = FT_2D(F,dir)
% Performs a 2D Fourier transform across the first two dimensions of F
%
% F: multidimensional array
% dir: direction; 'forward' or 'inverse'

if strcmp(dir,'forward')
    F = fftshift(fftshift(fft2(F),1),2);
end

if strcmp(dir,'inverse')
    F = ifft2(fftshift(fftshift(F,1),2));
end

end