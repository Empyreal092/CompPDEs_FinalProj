function [] = loglog_ordofconv(plot_input_ary,error_ary_mat,p,style)

arguments
    plot_input_ary
    error_ary_mat
    p
    style = ""
end


legend_nm = compose("%d",p);
loglog(plot_input_ary,min(error_ary_mat(2,:))/plot_input_ary(end)^p*plot_input_ary.^p,"k"+style,'DisplayName',"pow($\Delta$,"+legend_nm+")"); hold on

end

