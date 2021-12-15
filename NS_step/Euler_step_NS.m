function [omega_temp] = Euler_step_NS(omega_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend)
global L dt ext_sz finufft_interp

[u,v] = omega_2_uv(omega_temp);

x_depart = mod(x_mesh - u*dt , L);
y_depart = mod(y_mesh - v*dt , L);

if finufft_interp
    omega_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_depart,y_depart,"finufft");
else
    omega_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_depart,y_depart,"linear");
end


end

