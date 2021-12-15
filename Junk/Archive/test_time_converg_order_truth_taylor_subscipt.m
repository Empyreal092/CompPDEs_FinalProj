[~,~,~,omega_truth] = uvpomega_taylor(x_mesh,y_mesh,T);

%% 
omega_real = omega_final;

%%
error_norm_l1 = grid_func_norm(omega_real,omega_truth,"1");
error_norm_l2 = grid_func_norm(omega_real,omega_truth,"2");
error_norm_inf = grid_func_norm(omega_real,omega_truth,"inf");

error_ary_l1_omega = [error_ary_l1_omega error_norm_l1];
error_ary_l2_omega = [error_ary_l2_omega error_norm_l2];
error_ary_inf_omega = [error_ary_inf_omega error_norm_inf];

%%
% figure
% heatmap2d(omega_truth)
% pplot(8,0.8)
% figure
% heatmap2d(omega_real)
% pplot(8,0.8)
% 
% 1+1;