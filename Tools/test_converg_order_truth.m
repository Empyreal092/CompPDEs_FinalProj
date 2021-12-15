function [error_ary_mat] = test_converg_order_truth(data_final,data_truth,error_ary_mat)
% global Nx Ny L ext_sz

% x_ary = 0:L/Nx:L-L/Nx;
% y_ary = 0:L/Ny:L-L/Ny;
% [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
% x_ary_extend = 0-ext_sz*L/Nx:L/Nx:L-L/Nx+ext_sz*L/Nx;
% y_ary_extend = 0-ext_sz*L/Nx:L/Ny:L-L/Ny+ext_sz*L/Nx;
% [x_mesh_extend,y_mesh_extend] = meshgrid(x_ary_extend,y_ary_extend);

omega_truth = data_truth;
omega_real = data_final;

error_norm_l1 = grid_func_norm(omega_real,omega_truth,"1");
error_norm_l2 = grid_func_norm(omega_real,omega_truth,"2");
error_norm_inf = grid_func_norm(omega_real,omega_truth,"inf");

error_col = [error_norm_l1; error_norm_l2; error_norm_inf];
error_ary_mat = [error_ary_mat error_col];


end

