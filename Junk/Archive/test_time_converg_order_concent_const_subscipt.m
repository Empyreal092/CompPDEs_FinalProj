switch IC_type
    case "sinp"
        concent_truth = IC_sinp(x_mesh,y_mesh,p);
    case "3Gaussian"
        concent_truth = IC_3vort(x_mesh,y_mesh);
end
%% 
concent_real = concent_final;

error_norm_l1 = grid_func_norm(concent_real,concent_truth,"1");
error_norm_l2 = grid_func_norm(concent_real,concent_truth,"2");
error_norm_inf = grid_func_norm(concent_real,concent_truth,"inf");

error_ary_l1_c = [error_ary_l1_c error_norm_l1];
error_ary_l2_c = [error_ary_l2_c error_norm_l2];
error_ary_inf_c = [error_ary_inf_c error_norm_inf];