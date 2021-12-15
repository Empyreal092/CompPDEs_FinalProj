function [omega_temp,omega_temp_prev] = Trap_step_NS(omega_temp,omega_temp_prev,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend)
global L dt ext_sz finufft_interp

iter_count = 0;

[un,vn] = omega_2_uv(omega_temp);
if isnan(omega_temp_prev)
    [um,vm] = omega_2_uv(omega_temp);
else
    [um,vm] = omega_2_uv(omega_temp_prev);
end
up = 2*un-um; vp = 2*vn-vm;
un_iter = un; vn_iter = vn;
while iter_count < 1 % only one iteration
    x_depart_iter = mod(x_mesh - (un_iter+up)/2*dt,L);
    y_depart_iter = mod(y_mesh - (vn_iter+vp)/2*dt,L);
    % first order interpolation is enough here
    if finufft_interp
        un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_depart_iter,y_depart_iter,"finufft");
        vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_depart_iter,y_depart_iter,"finufft");
    else
        un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_depart_iter,y_depart_iter,"linear");
        vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_depart_iter,y_depart_iter,"linear");
    end

    
    iter_count = iter_count+1;
end
x_depart = mod((x_mesh - (un_iter+up)/2*dt),L);
y_depart = mod((y_mesh - (vn_iter+vp)/2*dt),L);

omega_temp_prev = omega_temp;
if finufft_interp
    omega_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_depart,y_depart,"finufft");
else
    omega_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_depart,y_depart,"cubic");
end


end

