function [data_final_prev,error_ary_mat] = test_converg_order_empirical(data_final,data_final_prev,error_ary_mat)
global Nx Ny L ext_sz

if ~isnan(data_final_prev)
    % Mx = 2^10+1; My = Mx;
    
%     omega_truth = omega_final_prev;
% %     omega_real = omega_final(1:2:end,1:2:end);
%     omega_real = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_final,x_mesh(1:2:end,1:2:end),y_mesh(1:2:end,1:2:end),"finufft");

    x_ary = 0:L/Nx:L-L/Nx;
    y_ary = 0:L/Ny:L-L/Ny;
    [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
    x_ary_extend = 0-ext_sz*L/Nx:L/Nx:L-L/Nx+ext_sz*L/Nx;
    y_ary_extend = 0-ext_sz*L/Nx:L/Ny:L-L/Ny+ext_sz*L/Nx;
    [x_mesh_extend,y_mesh_extend] = meshgrid(x_ary_extend,y_ary_extend);
    
    omega_truth = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,data_final_prev,x_mesh,y_mesh,"finufft");
    omega_real = data_final;
    
    error_norm_l1 = grid_func_norm(omega_real,omega_truth,"1");
    error_norm_l2 = grid_func_norm(omega_real,omega_truth,"2");
    error_norm_inf = grid_func_norm(omega_real,omega_truth,"inf");
    
%     error_ary_l1 = [error_ary_l1 error_norm_l1];
%     error_ary_l2 = [error_ary_l2 error_norm_l2];
%     error_ary_inf = [error_ary_inf error_norm_inf];

    error_col = [error_norm_l1; error_norm_l2; error_norm_inf];    
    error_ary_mat = [error_ary_mat error_col];
end

data_final_prev = data_final;

end

