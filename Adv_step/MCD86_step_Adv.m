function [tracer_temp,tracer_temp_prev]...
    = MCD86_step_Adv(tracer_temp,tracer_temp_prev,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,uh,vh)
global L dt ext_sz finufft_interp
iter_count = 0;

%% uh_hat
% uh = 1/8*(15*un-10*um+3*um2); vh = 1/8*(15*vn-10*vm+3*vm2);
uh_iter = uh; vh_iter = vh;

while iter_count < 2 % only 2 iteration
    x_h_iter = mod(x_mesh - uh_iter*dt/2,L);
    y_h_iter = mod(y_mesh - vh_iter*dt/2,L);
%     [uh_iter,vh_iter] = vel_func(x_h_iter,y_h_iter,T_curr+dt/2);
    if finufft_interp
        uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"finufft");
        vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"finufft");
    else
        if iter_count == 0
            uh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x_h_iter,y_h_iter,"linear");
            vh_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x_h_iter,y_h_iter,"linear");
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
    tracer_temp_hat = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_uhat,y_uhat,"finufft");
else
    tracer_temp_hat = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_uhat,y_uhat,"spline");
end

%% u_tilde
iter_count = 0;

un_iter = un; vn_iter = vn;

while iter_count < 2 % only 2 iteration
    x_n_iter = mod(x_mesh - un_iter*dt,L);
    y_n_iter = mod(y_mesh - vn_iter*dt,L);
%     [un_iter,vn_iter] = vel_func(x_n_iter,y_n_iter,T_curr);
    if finufft_interp
        un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"finufft");
        vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"finufft");
    else
        if iter_count == 0
            un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_n_iter,y_n_iter,"linear");
            vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_n_iter,y_n_iter,"linear");
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

if isnan(tracer_temp_prev)
    tracer_temp_tilde = tracer_temp_hat;
else
    if finufft_interp
        tracer_temp_tilde = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp_prev,x_utilde,y_utilde,"finufft");
    else
        tracer_temp_tilde = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp_prev,x_utilde,y_utilde,"spline");
    end
end

%% Bring them together
tracer_temp_prev = tracer_temp;
tracer_temp = 1/7*( 8*tracer_temp_hat-tracer_temp_tilde );
end

