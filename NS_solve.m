clear all
close all
addpath(genpath('Tools'))
addpath(genpath('NS_step'))
addpath(genpath('IC_n_Vel_Data'))

global L v0 Nx Ny dt ext_sz finufft_interp extrap_4th

%%
% IC_type = "Taylor"; N_resolve = 9;
IC_type = "3Vortices"; N_resolve = 33;

%%
% time_step_method = "Euler";
% time_step_method = "Trap";
% time_step_method = "MCD86"; % Fletcher p.218; McDonald, A. 1987
time_step_method = "RK4SL"; extrap_4th = false;
% time_step_method = "IF-RK4PS";

initialize_w_realdata = false;

finufft_interp = true;

disp("Time Step Method: "+time_step_method+"; Spectrual Interp: "+finufft_interp)

%%
if time_step_method == "IF-RK4PS"
    CFL_num = 1/8;
else
    CFL_num = 1.5;
end

%%
plot_timestep_numr = false;

if IC_type == "Taylor"
    if_test_converg_order_truth = true;
    if_test_converg_order_empiri = false;
else
    if_test_converg_order_truth = false;
    if_test_converg_order_empiri = true;
end

% some automatic options
fix_spacegrid_num = false;
if finufft_interp || time_step_method == "IF-RK4PS"
    fix_spacegrid_num = true;
end

%%
switch IC_type
    case "Taylor"
        L = 1;
        T = 0.25;
        v0 = 1;
    case "3Vortices"
        L = 2*pi;
        T = 1.5;
        v0 = 1;
end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    error_ary_mat = [];
    
    N_pow = [4:9];

    if finufft_interp
        N_pow = [1:6];
        N_ary = round(2.^N_pow);
    else
        N_ary = round(2.^N_pow)+1;
    end
    
    if if_test_converg_order_truth
        plot_input_ary = L./N_ary;
    else
        plot_input_ary = L./N_ary(1:end-1);
    end
else
    N_ary = 2^7;
end

if if_test_converg_order_empiri
    omega_final_prev = NaN;
end

if finufft_interp
    plot_input_ary = [];
end

for N = N_ary
    Nx = N; Ny = N;
    
    if finufft_interp
        Nt = N;
    else
        Nt = round( ((N/CFL_num)*(T/L)*v0)/2 )*2;
    end
    if finufft_interp && (N ~= N_ary(end)||IC_type == "Taylor")
        plot_input_ary = [plot_input_ary T/Nt];
    end
    
    if fix_spacegrid_num
        Nx = N_resolve; Ny = Nx;
    end
    dt = T/Nt;
    disp("Nx = "+Nx+"; Nt = "+Nt+"; CFL_Num: "+v0*Nx/Nt*(T/L));
    
    ext_sz = 5;
    x_ary = 0:L/Nx:L-L/Nx;
    y_ary = 0:L/Ny:L-L/Ny;
    [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
    x_ary_extend = 0-ext_sz*L/Nx:L/Nx:L-L/Nx+ext_sz*L/Nx;
    y_ary_extend = 0-ext_sz*L/Ny:L/Ny:L-L/Ny+ext_sz*L/Ny;
    [x_mesh_extend,y_mesh_extend] = meshgrid(x_ary_extend,y_ary_extend);
    
    switch IC_type
        case "Taylor"
            [~,~,IC_omega_real] = vel_taylor(x_mesh,y_mesh,0);
        case "3Vortices"
            IC_omega_real = IC_3vort(x_mesh,y_mesh);
    end
    
    %%
    T_curr = 0;
    omega_temp = IC_omega_real;
    
    if time_step_method == "Trap"
        omega_temp_prev = NaN;
    elseif time_step_method == "MCD86"
        omega_temp_prev = NaN;
        omega_temp_prev2 = NaN;
    elseif time_step_method == "RK4SL"
        if initialize_w_realdata
            [~,~,omega_temp_prev] = vel_taylor(x_mesh,y_mesh,-dt);
            [~,~,omega_temp_prev2] = vel_taylor(x_mesh,y_mesh,-2*dt);
            [~,~,omega_temp_prev3] = vel_taylor(x_mesh,y_mesh,-3*dt);
        else
            omega_temp_prev = NaN;
            omega_temp_prev2 = NaN;
            omega_temp_prev3 = NaN;
        end
    end
    
    while T_curr < T-dt/2
        
        switch time_step_method
            case "Euler"
                omega_temp = Euler_step_NS(omega_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
            case "Trap"
                [omega_temp,omega_temp_prev] = Trap_step_NS(omega_temp,omega_temp_prev,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
            case "MCD86"
                [omega_temp,omega_temp_prev,omega_temp_prev2] = ...
                    MCD86_step_NS(omega_temp,omega_temp_prev,omega_temp_prev2,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
            case "RK4SL"
                [omega_temp,omega_temp_prev,omega_temp_prev2,omega_temp_prev3] = ...
                    RK4SL_step_NS(omega_temp,omega_temp_prev,omega_temp_prev2,omega_temp_prev3,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend);
            case "IF-RK4PS"
                omega_temp = IF_RK4PS_step_NS(omega_temp);
        end
        
        T_curr = T_curr+dt;
        %%
        if plot_timestep_numr && N==N_ary(end)
            figure(99)
            heatmap2d(omega_temp,x_mesh,y_mesh); hold on
            if abs(T_curr-dt)<1e-10
                pplot(12,0.8)
            end
        end
    end
    if abs(T_curr-T)>dt/10
        disp("T_timestep not equal T!")
    end
    omega_final = omega_temp;
    
    %%
    if N==N_ary(end)
%         figure(101)
%         pplot(8,0.78,8)
%         heatmap2d(IC_omega_real,x_mesh,y_mesh); hold on
%         title("Initial Tracer Distrib $c(x,0)$")
%         xlabel("$x$"); ylabel("$y$")
%         pplot(8,0.78,8)
        
        figure(102)
        pplot(8,0.78,8)
        heatmap2d(omega_final,x_mesh,y_mesh); hold on
        title("$c(x,T)$; "+time_step_method)
        xlabel("$x$"); ylabel("$y$")
        pplot(8,0.78,8)
    end
    %%
    if if_test_converg_order_truth
        switch IC_type
            case "Taylor"
                [~,~,omega_truth] = vel_taylor(x_mesh,y_mesh,T);
                error_ary_mat = test_converg_order_truth(omega_final,omega_truth,error_ary_mat);
        end
    end
    if if_test_converg_order_empiri
        [omega_final_prev,error_ary_mat] = test_converg_order_empirical(omega_final,omega_final_prev,error_ary_mat);
    end
    
end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    figure(100)
    switch time_step_method
        case "Euler"
            loglog_ordofconv(plot_input_ary,error_ary_mat,1); hold on
        case "Trap"
            loglog_ordofconv(plot_input_ary,error_ary_mat,2); hold on
            loglog_ordofconv(plot_input_ary,error_ary_mat,1,"--");
        case "MCD86"
            loglog_ordofconv(plot_input_ary,error_ary_mat,3); hold on
            loglog_ordofconv(plot_input_ary,error_ary_mat,2,"--"); 
        case "RK4SL"
            loglog_ordofconv(plot_input_ary,error_ary_mat,4); hold on
            loglog_ordofconv(plot_input_ary,error_ary_mat,3,"--");
            loglog_ordofconv(plot_input_ary,error_ary_mat,2,":");
        case "IF-RK4PS"
            loglog_ordofconv(plot_input_ary,error_ary_mat,4); hold on
    end
    
    loglog(plot_input_ary,error_ary_mat(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
    loglog(plot_input_ary,error_ary_mat(2,:),'b^','DisplayName','$c,\ell^2$')
    loglog(plot_input_ary,error_ary_mat(3,:),'bs','DisplayName','$c$, uniform')
    
    xlim([plot_input_ary(end) plot_input_ary(1)])
    ylabel('error'), 
    if finufft_interp
        title(["IC: "+IC_type+"; Method: "+time_step_method+" w/ FINUFFT"])
        xlabel('$\Delta t$')
    else
        title(["IC: "+IC_type+"; Method: "+time_step_method+"; CFL $\approx $"+CFL_num])
        xlabel('$\Delta x$')
    end
    
    pplot(8,0.8,8)
    legend('Location','best','NumColumns',1)
    hold off
end

%%
if finufft_interp
    nm_tag = "_finu_";
else
    nm_tag = "_";
end

figure(100)
savefig("latex/figs/"+"nonlin_conv_order_"+IC_type+nm_tag+time_step_method)
figure(102)
savefig("latex/figs/"+"nonlin_omega_final_"+IC_type+nm_tag+time_step_method)

