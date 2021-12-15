% Mx = 2^10+1; My = Mx;
x_dense = 0:L/Mx:L-L/Mx; y_dense = 0:L/My:L-L/My;
[x_mesh_d,y_mesh_d] = meshgrid(x_dense,y_dense);

[u_truth,~,~,~] = taylor_exact(x_mesh_d,y_mesh_d,T);
%% 
concent_final_interpn = zero_pad(concent_final,Mx,My);
u_real = real(ifft2_n(concent_final_interpn));

error_norm_l1 = grid_func_norm(u_real,u_truth,"1");
error_norm_l2 = grid_func_norm(u_real,u_truth,"2");
error_norm_inf = grid_func_norm(u_real,u_truth,"inf");

error_ary_l1_uv = [error_ary_l1_uv error_norm_l1];
error_ary_l2_uv = [error_ary_l2_uv error_norm_l2];
error_ary_inf_uv = [error_ary_inf_uv error_norm_inf];

% error_norm_l1 = grid_func_norm(omega_timestep_real,omega_truth,"1");
% error_norm_l2 = grid_func_norm(omega_timestep_real,omega_truth,"2");
% error_norm_inf = grid_func_norm(omega_timestep_real,omega_truth,"inf");
% 
% error_ary_l1_omega = [error_ary_l1_omega error_norm_l1];
% error_ary_l2_omega = [error_ary_l2_omega error_norm_l2];
% error_ary_inf_omega = [error_ary_inf_omega error_norm_inf];