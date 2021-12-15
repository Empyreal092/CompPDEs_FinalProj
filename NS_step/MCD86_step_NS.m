function [omega_temp,omega_temp_prev,omega_temp_prev2] = MCD86_step_NS(omega_temp,omega_temp_prev,omega_temp_prev2,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend)
global L dt ext_sz finufft_interp

iter_count = 0;

[un,vn] = omega_2_uv(omega_temp);

if isnan(omega_temp_prev)
    omega_temp_prev = omega_temp;
    dt = dt/16;
    for i = 1:16
        omega_temp = IF_RK4PS_step_NS(omega_temp);
    end
    dt = dt*16;
elseif isnan(omega_temp_prev2)
    omega_temp_prev2 = omega_temp_prev;
    [omega_temp,omega_temp_prev] = Trap_step_NS(omega_temp,omega_temp_prev,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
else
    [um,vm] = omega_2_uv(omega_temp_prev);
    [um2,vm2] = omega_2_uv(omega_temp_prev2);
    %% MCD86 step
    %% uh_hat
    uh = 1/8*(15*un-10*um+3*um2); vh = 1/8*(15*vn-10*vm+3*vm2);
    uh_iter = uh; vh_iter = vh;
    
    while iter_count < 2 % only 2 iteration
        x_h_iter = mod(x_mesh - uh_iter*dt/2,L);
        y_h_iter = mod(y_mesh - vh_iter*dt/2,L);
        if iter_count == 0
            if finufft_interp
                uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"finufft");
                vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"finufft");
            else
                uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"linear");
                vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"linear");
            end
            
        else
            if finufft_interp
                uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"finufft");
                vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"finufft");
            else
                uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"cubic");
                vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"cubic");
            end
            
        end
        iter_count = iter_count+1;
    end
    uh_hat = uh_iter; vh_hat = vh_iter;
    
    x_uhat = mod((x_mesh - uh_hat*dt),L);
    y_uhat = mod((y_mesh - vh_hat*dt),L);
    
    if finufft_interp
        omega_temp_hat = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_uhat,y_uhat,"finufft");
    else
        omega_temp_hat = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp,x_uhat,y_uhat,"spline");
    end
    
    
    %% u_tilde
    iter_count = 0;
    
    un_iter = un; vn_iter = vn;
    
    while iter_count < 2 % only 2 iteration
        x_n_iter = mod(x_mesh - un_iter*dt,L);
        y_n_iter = mod(y_mesh - vn_iter*dt,L);
        if iter_count == 0
            if finufft_interp
                un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"finufft");
                vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"finufft");
            else
                un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"linear");
                vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"linear");
            end
            
        else
            if finufft_interp
                un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"finufft");
                vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"finufft");
            else
                un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"cubic");
                vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"cubic");
            end
            
        end
        iter_count = iter_count+1;
    end
    u_tilde = un_iter; v_tilde = vn_iter;
    
    x_utilde = mod((x_mesh - u_tilde*dt*2),L);
    y_utilde = mod((y_mesh - v_tilde*dt*2),L);
    
    if finufft_interp
        omega_temp_tilde = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp_prev,x_utilde,y_utilde,"finufft");
    else
        omega_temp_tilde = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_temp_prev,x_utilde,y_utilde,"spline");
    end
    
    %% Bring them together
    omega_temp_prev2 = omega_temp_prev;
    omega_temp_prev = omega_temp;
    omega_temp = 1/7*( 8*omega_temp_hat-omega_temp_tilde );
end



end

