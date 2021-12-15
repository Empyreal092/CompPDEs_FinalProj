% Mx = 2^10+1; My = Mx;
x_dense = 0:L/Mx:L-L/Mx; y_dense = 0:L/My:L-L/My;
[x_mesh_d,y_mesh_d] = meshgrid(x_dense,y_dense);

[u_truth,v_truth,~,omega_truth] = taylor_exact(x_mesh_d,y_mesh_d,T);
%% 
omega_final_interpn = zero_pad(omega_final,Mx,My);
k_dense = [0:(Mx-1)/2 -(Mx-1)/2:-1]*2*pi/L; l_dense = [0:(My-1)/2 -(My-1)/2:-1]*2*pi/L;
[k_mesh_d,l_mesh_d] = meshgrid(k_dense,l_dense);
inv_laplace_mat = inv_laplace_make(k_dense,l_dense);
psi = inv_laplace_mat.*omega_final_interpn;
u_four = -1i*psi.*l_mesh_d;
v_four =  1i*psi.*k_mesh_d;
u_real = v0+real(ifft2_n(u_four)); v_real = v0+real(ifft2_n(v_four));
omega_timestep_real = real(ifft2_n(omega_final_interpn));

error_norm_l1 = grid_func_norm([grid_func_norm(u_real,u_truth,"1") grid_func_norm(v_real,v_truth,"1")],zeros(2,1),"1");
error_norm_l2 = grid_func_norm([grid_func_norm(u_real,u_truth,"2") grid_func_norm(v_real,v_truth,"2")],zeros(2,1),"2");
error_norm_inf = grid_func_norm([grid_func_norm(u_real,u_truth,"inf") grid_func_norm(v_real,v_truth,"inf")],zeros(2,1),"inf");

error_ary_l1_uv = [error_ary_l1_uv error_norm_l1];
error_ary_l2_uv = [error_ary_l2_uv error_norm_l2];
error_ary_inf_uv = [error_ary_inf_uv error_norm_inf];

error_norm_l1 = grid_func_norm(omega_timestep_real,omega_truth,"1");
error_norm_l2 = grid_func_norm(omega_timestep_real,omega_truth,"2");
error_norm_inf = grid_func_norm(omega_timestep_real,omega_truth,"inf");

error_ary_l1_omega = [error_ary_l1_omega error_norm_l1];
error_ary_l2_omega = [error_ary_l2_omega error_norm_l2];
error_ary_inf_omega = [error_ary_inf_omega error_norm_inf];