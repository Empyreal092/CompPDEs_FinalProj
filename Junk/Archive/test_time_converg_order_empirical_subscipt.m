if ~isnan(concent_final_prev)
    % Mx = 2^10+1; My = Mx;
    concent_truth = concent_final_prev;
    %%
    concent_real = concent_final(1:2:end,1:2:end);
    
    error_norm_l1 = grid_func_norm(concent_real,concent_truth,"1");
    error_norm_l2 = grid_func_norm(concent_real,concent_truth,"2");
    error_norm_inf = grid_func_norm(concent_real,concent_truth,"inf");
    
    error_ary_l1_c = [error_ary_l1_c error_norm_l1];
    error_ary_l2_c = [error_ary_l2_c error_norm_l2];
    error_ary_inf_c = [error_ary_inf_c error_norm_inf];
end

%% 
concent_final_prev = concent_final;