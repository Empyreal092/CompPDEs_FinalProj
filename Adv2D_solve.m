clear all
close all
addpath(genpath('Tools'))
addpath(genpath('AdvStep'))
addpath(genpath('IC_n_Vel_Data'))

global L v0 Nx Ny dt ext_sz finufft_interp

%%
% vel_type = "Constant";
vel_type = "Taylor";

% IC_type = "sinp"; p = 2; N_resolve = 9;
IC_type = "3Gaussian"; N_resolve = 81;

%%
% time_step_method = "Euler";
% time_step_method = "Trap";
% time_step_method = "MCD86"; % Fletcher p.218; McDonald, A. 1987
time_step_method = "RK4SL";
% time_step_method = "IF-RK4PS";

finufft_interp = true;

disp("Time Step Method: "+time_step_method+"; Spectrual Interp: "+finufft_interp)

%% 
switch vel_type
    case "Constant"
        if time_step_method == "IF-RK4PS"
            CFL_num = 1/5;
        else
            CFL_num = 8;
        end
    case "Taylor"
        CFL_num = 1.5;
        if time_step_method == "IF-RK4PS"
            CFL_num = 0.5;
        end
end
%%
plot_timestep_numr = false;

if_test_converg_order_truth = false;
if_test_converg_order_empiri = true;

% some automatic options
fix_timesteps_num = false;
if vel_type == "Constant"
    fix_timesteps_num = true;
end
fix_spacegrid_num = false;
if finufft_interp || time_step_method == "IF-RK4PS"
    fix_spacegrid_num = true;
end

%%
switch IC_type
    case "sinp"
        L = 1;
    case "3Gaussian"
        L = 2*pi;
end

switch vel_type
    case "Constant"
        T = 2;
        v0 = L;
        vel_func = @vel_constant;
    case "Taylor"
        T = 0.5;
        v0 = 1;
        vel_func = @vel_taylor;
end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    error_ary_mat = [];
    
    N_pow = [5:10];
    N_ary = round(2.^N_pow)+1;
    
    if if_test_converg_order_truth
        plot_input_ary = L./N_ary;
    else
        plot_input_ary = L./N_ary(1:end-1);
    end
else
    N_ary = 2^7;
end

if if_test_converg_order_empiri
    tracer_final_prev = NaN;
end

for N = N_ary
    Nx = N; Ny = N;
    Nt = round( ((N/CFL_num)*(T/L)*v0)/2 )*2;

    if fix_timesteps_num
        Nt = 15;
    end
    if fix_spacegrid_num
        Nx = N_resolve; Ny = Nx;
    end
    dt = T/Nt;
    disp("Nx = "+Nx+"; Nt = "+Nt+"; CFL_Num: "+v0*Nx/Nt*(T/L));
    
    ext_sz = 5;
    x_ary = 0:L/Nx:L-L/Nx; y_ary = 0:L/Ny:L-L/Ny;
    [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
    x_ary_extend = 0-ext_sz*L/Nx:L/Nx:L-L/Nx+ext_sz*L/Nx; y_ary_extend = 0-ext_sz*L/Nx:L/Ny:L-L/Ny+ext_sz*L/Nx;
    [x_mesh_extend,y_mesh_extend] = meshgrid(x_ary_extend,y_ary_extend);
    
    switch IC_type
        case "sinp"
            IC_tracer_real = IC_sinp(x_mesh,y_mesh,p);
        case "3Gaussian"
            IC_tracer_real = IC_3vort(x_mesh,y_mesh);
    end
    %%
    T_curr = 0;
    tracer_temp = IC_tracer_real;
    
    if time_step_method == "MCD86"
        tracer_temp_prev = NaN;
    end
    
    while T_curr < T-dt/2
        
        switch time_step_method
            case "Euler"
                [un,vn] = vel_func(x_mesh,y_mesh,T_curr);
                tracer_temp = Euler_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn);
            case "Trap"
                [un,vn] = vel_func(x_mesh,y_mesh,T_curr);
                [um,vm] = vel_func(x_mesh,y_mesh,T_curr-dt);
                tracer_temp = Trap_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,um,vm);
             case "MCD86"
                 [un,vn] = vel_func(x_mesh,y_mesh,T_curr); 
                 [uh,vh] = vel_func(x_mesh,y_mesh,T_curr+dt/2);
%                  [um2,vm2] = vel_func(x_mesh,y_mesh,T_curr-2*dt);
                 [tracer_temp,tracer_temp_prev]...
                     = MCD86_step_Adv(tracer_temp,tracer_temp_prev,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,uh,vh);
            case "RK4SL"
                 [un,vn] = vel_func(x_mesh,y_mesh,T_curr); 
                 [uh,vh] = vel_func(x_mesh,y_mesh,T_curr+dt/2);
                 [up,vp] = vel_func(x_mesh,y_mesh,T_curr+dt);
                 tracer_temp = RK4SL_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,uh,vh,up,vp);
            case "IF-RK4PS"
                 [~,~,omega_n] = vel_func(x_mesh,y_mesh,T_curr); 
                 [~,~,omega_h] = vel_func(x_mesh,y_mesh,T_curr+dt/2);
                 [~,~,omega_p] = vel_func(x_mesh,y_mesh,T_curr+dt);
                 tracer_temp = IF_RK4PS_step_Adv(tracer_temp,omega_n,omega_h,omega_p);
        end
        
        T_curr = T_curr+dt;
        %%
        if plot_timestep_numr && N==N_ary(end)
            figure(99)
            heatmap2d(tracer_temp,x_mesh,y_mesh); hold on
            if abs(T_curr-dt)<1e-10
                pplot(12,0.8)
            end
        end
    end
    if abs(T_curr-T)>dt/10
        disp("T_timestep not equal T!")
    end
    tracer_final = tracer_temp;
    
    %% 
    if N==N_ary(end)
        figure(101)
        pplot(8,0.78,8)
        heatmap2d(IC_tracer_real,x_mesh,y_mesh); hold on
        title("Initial Tracer Distrib $c(x,0)$")
        xlabel("$x$"); ylabel("$y$")
        pplot(8,0.78,8)
        
        figure(102)
        pplot(8,0.78,8)
        heatmap2d(tracer_final,x_mesh,y_mesh); hold on
        title("$c(x,T)$; "+time_step_method)
        xlabel("$x$"); ylabel("$y$")
        pplot(8,0.78,8)
    end
    %%
    if if_test_converg_order_truth && vel_type == "Constant"
        tracer_truth = IC_tracer_real;
        error_ary_mat = test_converg_order_truth(tracer_final,tracer_truth,error_ary_mat);
    end
    if if_test_converg_order_empiri
        [tracer_final_prev,error_ary_mat] = test_converg_order_empirical(tracer_final,tracer_final_prev,error_ary_mat);
    end

end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    figure(100)
    
    if vel_type == "Constant" && time_step_method ~= "IF-RK4PS" && finufft_interp
        semilogy(plot_input_ary,error_ary_mat(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
        semilogy(plot_input_ary,error_ary_mat(2,:),'b^','DisplayName','$c,\ell^2$')
        semilogy(plot_input_ary,error_ary_mat(3,:),'bs','DisplayName','$c$, uniform')
    else
        loglog(plot_input_ary,error_ary_mat(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
        loglog(plot_input_ary,error_ary_mat(2,:),'b^','DisplayName','$c,\ell^2$')
        loglog(plot_input_ary,error_ary_mat(3,:),'bs','DisplayName','$c$, uniform')
        
        switch time_step_method
            case "Euler"
                loglog_ordofconv(plot_input_ary,error_ary_mat,1)
            case "Trap"
                loglog_ordofconv(plot_input_ary,error_ary_mat,2)
            case "MCD86"
                loglog_ordofconv(plot_input_ary,error_ary_mat,3)
            case {"RK4SL","IF-RK4PS"}
                loglog_ordofconv(plot_input_ary,error_ary_mat,4)
        end
    end

    xlim([plot_input_ary(end) plot_input_ary(1)])
    
    ylabel('error'), xlabel('$\Delta x$')
    switch vel_type
        case "Constant"
%             title([vel_type+" Vel; "+"IC: "+IC_type])
    end

    pplot(8,0.8,8)
    legend('Location','best','NumColumns',1)
    title("Method: "+time_step_method)
    hold off
end

%%
if finufft_interp
    figure(100)
    savefig("latex/figs/"+"conv_order_finu_"+time_step_method)
%     figure(101)
%     savefig("latex/figs/"+"c_init_finu_"+time_step_method)
%     figure(102)
%     savefig("latex/figs/"+"c_final_finu_"+time_step_method)
else
%     figure(100)
%     savefig("latex/figs/"+"conv_order_"+time_step_method)
%     figure(101)
%     savefig("latex/figs/"+"c_init_"+time_step_method)
%     figure(102)
%     savefig("latex/figs/"+"c_final_"+time_step_method)
end
