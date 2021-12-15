if ~isnan(omega_final_prev)
    % Mx = 2^10+1; My = Mx;
    
%     omega_truth = omega_final_prev;
% %     omega_real = omega_final(1:2:end,1:2:end);
%     omega_real = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_final,x_mesh(1:2:end,1:2:end),y_mesh(1:2:end,1:2:end),"finufft");
    
    omega_truth = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,omega_final_prev,x_mesh,y_mesh,"finufft");
    omega_real = omega_final;
    
    error_norm_l1 = grid_func_norm(omega_real,omega_truth,"1");
    error_norm_l2 = grid_func_norm(omega_real,omega_truth,"2");
    error_norm_inf = grid_func_norm(omega_real,omega_truth,"inf");
    
    error_ary_l1_omega = [error_ary_l1_omega error_norm_l1];
    error_ary_l2_omega = [error_ary_l2_omega error_norm_l2];
    error_ary_inf_omega = [error_ary_inf_omega error_norm_inf];
end

omega_final_prev = omega_final;