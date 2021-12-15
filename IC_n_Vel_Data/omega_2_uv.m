function [u,v] = omega_2_uv(omega)
global L v0

omega_four = fft2_n(omega);

[Ny,Nx] = size(omega_four);
k_ary = [0:(Nx-1)/2 -(Nx-1)/2:-1]*2*pi/L;
l_ary = [0:(Ny-1)/2 -(Ny-1)/2:-1]*2*pi/L;
[k_mesh,l_mesh] = meshgrid(k_ary,l_ary);

inv_laplace_mat = inv_laplace_make(k_ary,l_ary);

psi = inv_laplace_mat.*omega_four;
u_four = -1i*psi.*l_mesh;
v_four =  1i*psi.*k_mesh;

u = ifft2_n(u_four)+v0;
v = ifft2_n(v_four)+v0;

end

