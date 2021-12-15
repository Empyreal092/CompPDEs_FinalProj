function [omega_temp,omega_temp_prev,omega_temp_prev2,omega_temp_prev3] = ...
    RK4SL_step_NS(omega_temp,omega_temp_prev,omega_temp_prev2,omega_temp_prev3,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend)
global L dt ext_sz extrap_4th finufft_interp


if isnan(omega_temp_prev3)
    omega_temp_prev3 = omega_temp_prev2;
    [omega_temp,omega_temp_prev,omega_temp_prev2] = ...
        MCD86_step_NS(omega_temp,omega_temp_prev,omega_temp_prev2,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
else
    %% Extrapolate Velocity
    [un,vn] = omega_2_uv(omega_temp);
    [um,vm] = omega_2_uv(omega_temp_prev);
    [um2,vm2] = omega_2_uv(omega_temp_prev2);
    [um3,vm3] = omega_2_uv(omega_temp_prev3);
    
    if extrap_4th
        uh = 1/16*(35*un-35*um+21*um2-5*um3); vh = 1/16*(35*vn-35*vm+21*vm2-5*vm3);
        up = (4*un-6*um+4*um2-um3);     vp = (4*vn-6*vm+4*vm2-vm3);
    else
        uh = 1/8*(15*un-10*um+3*um2); vh = 1/8*(15*vn-10*vm+3*vm2);
        up = (3*un-3*um+um2);     vp = (3*vn-3*vm+vm2);
    end

    %%
    if finufft_interp
        in_method_interp = "finufft";
    else
        in_method_interp = "spline";
    end

    x1 = mod(x_mesh,L); y1 = mod(y_mesh,L);
    % [u1,v1,~,~] = uvpomega_taylor(x1,y1,T_timestep+dt);
    u1 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,up,x1,y1,in_method_interp);
    v1 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vp,x1,y1,in_method_interp);
    
    x2 = mod(x_mesh-dt/2*u1,L); y2 = mod(y_mesh-dt/2*v1,L);
    % [u2,v2,~,~] = uvpomega_taylor(x2,y2,T_timestep+dt/2);
    u2 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x2,y2,in_method_interp);
    v2 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x2,y2,in_method_interp);
    
    x3 = mod(x_mesh-dt/2*u2,L); y3 = mod(y_mesh-dt/2*v2,L);
    % [u3,v3,~,~] = uvpomega_taylor(x3,y3,T_timestep+dt/2);
    u3 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x3,y3,in_method_interp);
    v3 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x3,y3,in_method_interp);
    
    x4 = mod(x_mesh-dt*u3,L); y4 = mod(y_mesh-dt*v3,L);
    % [u4,v4,~,~] = uvpomega_taylor(x4,y4,T_timestep);
    u4 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x4,y4,in_method_interp);
    v4 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x4,y4,in_method_interp);
    
    um = (u1 + 2*u2 + 2*u3 + u4)/6;  vm = (v1 + 2*v2 + 2*v3 + v4)/6;
    
    x_depart = mod((x_mesh - um*dt),L);
    y_depart = mod((y_mesh - vm*dt),L);
    
    omega_temp_prev2 = omega_temp_prev;
    omega_temp_prev = omega_temp;
    omega_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_depart,y_depart,"finufft");
end


end

