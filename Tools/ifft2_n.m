%standard 2D fourier inverse transformation
function res=ifft2_n(x)
[Nx,Ny]=size(x);
res=ifft2(x)*(Nx*Ny); % correct normalization
end