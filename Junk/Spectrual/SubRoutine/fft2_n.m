%standard 2D fourier transformation
function res=fft2_n(x)
[Nx,Ny]=size(x);
res=fft2(x)/(Nx*Ny); % correct normalization
end