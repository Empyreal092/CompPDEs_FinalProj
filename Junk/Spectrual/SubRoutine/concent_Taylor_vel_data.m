function [omega_data_four,pex_pres_four] = concent_Taylor_vel_data(x_mesh,y_mesh,T_timestep)
global L

[Ny,Nx] = size(x_mesh);
k_ary = [0:(Nx-1)/2 -(Nx-1)/2:-1]*2*pi/L;
l_ary = [0:(Ny-1)/2 -(Ny-1)/2:-1]*2*pi/L;
[k_mesh,~] = meshgrid(k_ary,l_ary);

[~,~,p_real,omega_real] = taylor_exact(x_mesh,y_mesh,T_timestep);
omega_data_four = fft2_n(omega_real);

p_four = fft2_n(p_real);
pex_pres_four = 1i*p_four.*k_mesh;

end

