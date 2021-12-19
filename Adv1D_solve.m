clear all
close all
addpath(genpath('Tools'))
% addpath(genpath('AdvStep'))
addpath(genpath('IC_n_Vel_Data'))

% global L v0 Nx Ny dt ext_sz finufft_interp
global L

%%
vel_type = "Constant";

% IC_type = "sinp"; p = 40;
IC_type = "step";

%%
time_step_method = "Euler";

% interp_method = "linear";
% interp_method = "cubic";
% interp_method = "pchip";
% interp_method = "makima";
% interp_method = "spline";
interp_method = "finufft";

disp("Interpolation Method: "+interp_method)

%% 
%%
plot_timestep_numr = true;
if_test_converg_order_truth = true;

% some automatic options
fix_timesteps_num = true;
fix_spacegrid_num = false;

%%
L = 1;
T = 1;
v0 = L;

%%
switch IC_type
    case "sinp"
        N_plot = 2^4;
    case "step"
        N_plot = 2^6;
end
    
if if_test_converg_order_truth
    error_ary_mat = [];
    
    N_pow = [3:10];
    N_ary = round(2.^N_pow);
    
    plot_input_ary = L./N_ary;
else
    N_ary = N_plot;
end

for N = N_ary
    Nx = N;
    Nt = 7;
    
    dt = T/Nt;
    disp("Nx = "+Nx+"; Nt = "+Nt+"; CFL_Num: "+v0*Nx/Nt*(T/L));
    
    ext_sz = 5;
    x_ary = 0:L/Nx:L-L/Nx; 
    x_mesh = x_ary;
    x_ary_extend = 0-ext_sz*L/Nx:L/Nx:L-L/Nx+ext_sz*L/Nx;
    x_mesh_extend = x_ary_extend;
    
    switch IC_type
        case "sinp"
            IC_tracer_real = IC_sinp1D(x_mesh,p);
        case "step"
            IC_tracer_real = IC_step(x_mesh);
    end
    %%
    T_curr = 0;
    tracer_temp = IC_tracer_real;
    
    if plot_timestep_numr && N==N_plot
        figure(99)
        x_IC_real_plot = 0:L/128:L-L/128;
        switch IC_type
            case "sinp"
                IC_real_plot = IC_sinp1D(x_IC_real_plot,p);
            case "step"
                IC_real_plot = IC_step(x_IC_real_plot);
        end
        plot([x_IC_real_plot x_IC_real_plot(1)+L],[IC_real_plot IC_real_plot(1)],'k-','DisplayName',"Truth"); hold on
        pplot(8,0.8,8)
    end

    while T_curr < T-dt/2
        un = v0;
        x_depart = mod(x_mesh - un*dt , L);
        tracer_temp = interp1_more(ext_sz,x_mesh_extend,tracer_temp,x_depart,interp_method);
        
        T_curr = T_curr+dt;
        %%
        if plot_timestep_numr
            figure(99)
            if T_curr > T-dt/2
                if N==2^4
                    tracer_temp_course = tracer_temp;
                    x_mesh_coarse = x_mesh;
                elseif N==2^6
                    plot([x_mesh x_mesh(1)+L],[tracer_temp tracer_temp(1)],'r-*','DisplayName',"$Nx=$"+64); hold on
                    plot([x_mesh_coarse x_mesh_coarse(1)+L],[tracer_temp_course tracer_temp_course(1)],'b-*','DisplayName',"$Nx=$"+16); hold on
                end
            end
        end
    end
    if abs(T_curr-T)>dt/10
        disp("T_timestep not equal T!")
    end
    tracer_final = tracer_temp;
    
    %%
    if if_test_converg_order_truth && vel_type == "Constant"
        tracer_truth = IC_tracer_real;
        error_ary_mat = test_converg_order_truth(tracer_final,tracer_truth,error_ary_mat);
    end

end

%%
figure(99)
ylim([-0.2 1.2])
ylabel('$c$'), xlabel('$x$')
title("$Nt=$"+Nt+"; Interp Method: "+interp_method)

legend('Location','northeast','NumColumns',1)
hold off

%%
if if_test_converg_order_truth
    figure(100)
    
    if interp_method == "finufft" && IC_type ~= "step"
        semilogy(plot_input_ary,error_ary_mat(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
        semilogy(plot_input_ary,error_ary_mat(2,:),'b^','DisplayName','$c,\ell^2$')
        semilogy(plot_input_ary,error_ary_mat(3,:),'bs','DisplayName','$c$, uniform')
    else
        loglog(plot_input_ary,error_ary_mat(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
        loglog(plot_input_ary,error_ary_mat(2,:),'b^','DisplayName','$c,\ell^2$')
        loglog(plot_input_ary,error_ary_mat(3,:),'bs','DisplayName','$c$, uniform')
        
        switch interp_method
            case "linear"
                loglog_ordofconv(plot_input_ary,error_ary_mat,2)
            case {"cubic"}
                loglog_ordofconv(plot_input_ary,error_ary_mat,3)
            case {"pchip","makima"}
                loglog_ordofconv(plot_input_ary,error_ary_mat(2:3,:),2,"--")
                loglog_ordofconv(plot_input_ary,error_ary_mat([2 1],:),3)
            case "spline"
                loglog_ordofconv(plot_input_ary,error_ary_mat,4)
            case "finufft"
                loglog_ordofconv(plot_input_ary,error_ary_mat,1)
        end
    end
        
    xlim([plot_input_ary(end) plot_input_ary(1)])
    
    ylabel('error'), xlabel('$\Delta x$')
    title("$Nt=$"+Nt+"; Interp Method: "+interp_method)

    pplot(8,0.8,8)
    legend('Location','best','NumColumns',1)
    hold off
end

%%
% file_nm = "1D_cons_"+interp_method;
% file_nm = "1D_step_cons_"+interp_method;
% figure(99)
% savefig("latex/figs/"+file_nm+"_sol")
% figure(100)
% savefig("latex/figs/"+file_nm+"_convord")
