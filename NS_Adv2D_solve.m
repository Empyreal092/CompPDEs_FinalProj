clear all
close all
addpath(genpath('Tools'))
addpath(genpath('NS_step'))
addpath(genpath('Adv_step'))
addpath(genpath('IC_n_Vel_Data'))

global L v0 Nx Ny dt ext_sz finufft_interp

%%
% vel_type = "Taylor";
% IC_type = "3Gaussian"; N_resolve = 81;

vel_type = "3Vortices";
IC_type = "sinp"; p = 2; N_resolve = 49;

%%
NS_time_step_method = "IF-RK4PS";
adv_time_step_method = "RK4SL";

finufft_interp = true;

disp("Time Step Method for Advection: "+adv_time_step_method+"; for NS: "+NS_time_step_method+"; Spectrual Interp: "+finufft_interp)

%%
% NS_CFL_num = 1/8;
adv_CFL_num = 2;

%%
plot_timestep_numr = false;

% some automatic options
if vel_type == "Taylor" && IC_type == "3Gaussian"
    if_test_converg_order_truth = true;
    if_test_converg_order_empiri = false;
else
    if_test_converg_order_truth = false;
    if_test_converg_order_empiri = true;
end

fix_timesteps_num = false;
    fix_spacegrid_num = true;
%%
switch IC_type
    case "sinp"
        L = 2*pi;
    case "3Gaussian"
        L = 2*pi;
end

switch vel_type
    case "3Vortices"
        T = 4;
        v0 = 1;
    case "Taylor"
        T = 0.5;
        v0 = 1;
end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    error_ary_mat_trac = [];
    error_ary_mat_vel =  [];
    
%     N_pow = [1:5];
    N_pow = [4:8];
    N_ary = round(2.^N_pow);
    
    if if_test_converg_order_truth
        plot_input_ary = T./N_ary;
        load('tracer_taylor_3vort_truth.mat')
    else
        plot_input_ary = T./N_ary(1:end-1);
    end
else
    N_ary = 2^7;
end

if if_test_converg_order_empiri
    tracer_final_prev = NaN;
    omega_final_prev = NaN;
end

for N = N_ary
    if fix_spacegrid_num
        Nx = N_resolve; Ny = Nx;
    end
%     Nt = round( ((Nx/adv_CFL_num)*(T/L)*v0*N)/2 )*2;
    Nt = N;
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
    switch vel_type
        case "Taylor"
            [~,~,IC_omega_real] = vel_taylor(x_mesh,y_mesh,0);
        case "3Vortices"
            IC_omega_real = IC_3vort(x_mesh,y_mesh);
    end
    %%
    T_curr = 0;
    tracer_temp = IC_tracer_real;
    omega_temp = IC_omega_real;
    
    while T_curr < T-dt/2
        
        switch adv_time_step_method
            case "RK4SL"
                [un,vn] = omega_2_uv(omega_temp);
                
                NS_to_adv_ratio = 2^2;
                dt = dt/NS_to_adv_ratio;
                for i = 1:NS_to_adv_ratio
                    omega_temp = IF_RK4PS_step_NS(omega_temp);
                    if i == NS_to_adv_ratio/2
                        [uh,vh] = omega_2_uv(omega_temp);
                    end
                    if i == NS_to_adv_ratio
                        [up,vp] = omega_2_uv(omega_temp);
                    end
                end
                dt = dt*NS_to_adv_ratio;
                
                tracer_temp = RK4SL_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,uh,vh,up,vp);
        end
        
        T_curr = T_curr+dt;
        %%
        if plot_timestep_numr && N==N_ary(end)
            figure(98)
            heatmap2d(tracer_temp,x_mesh,y_mesh); hold on
            if abs(T_curr-dt)<1e-10
                pplot(12,0.8)
            end
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
    tracer_final = tracer_temp;
    
    %%
    if if_test_converg_order_truth && vel_type == "Taylor" && IC_type == "3Gaussian"
        tracer_truth = tracer_taylor_3vort_truth;
        [~,~,omega_truth] = vel_taylor(x_mesh,y_mesh,T);
        error_ary_mat_trac = test_converg_order_truth(tracer_final,tracer_truth,error_ary_mat_trac);
        error_ary_mat_vel = test_converg_order_truth(omega_final,omega_truth,error_ary_mat_vel);
    end
    if if_test_converg_order_empiri
        [tracer_final_prev,error_ary_mat_trac] = test_converg_order_empirical(tracer_final,tracer_final_prev,error_ary_mat_trac);
        [omega_final_prev,error_ary_mat_vel] = test_converg_order_empirical(omega_final,omega_final_prev,error_ary_mat_vel);
    end
    
end

%%
if if_test_converg_order_truth || if_test_converg_order_empiri
    figure(100)
    
    switch adv_time_step_method
        case "Euler"
            loglog_ordofconv(plot_input_ary,error_ary_mat_trac,1)
        case "Trap"
            loglog_ordofconv(plot_input_ary,error_ary_mat_trac,2)
        case "MCD86"
            loglog_ordofconv(plot_input_ary,error_ary_mat_trac,3)
        case {"RK4SL","IF-RK4PS"}
            loglog_ordofconv(plot_input_ary,error_ary_mat_trac,4)
    end
    
    loglog(plot_input_ary,error_ary_mat_trac(1,:),'bo','DisplayName','$c,\ell^1$'); hold on
    loglog(plot_input_ary,error_ary_mat_trac(2,:),'b^','DisplayName','$c,\ell^2$')
    loglog(plot_input_ary,error_ary_mat_trac(3,:),'bs','DisplayName','$c$, uniform')
    loglog(plot_input_ary,error_ary_mat_vel(1,:),'ro','DisplayName','$\omega,\ell^1$')
    loglog(plot_input_ary,error_ary_mat_vel(2,:),'r^','DisplayName','$\omega,\ell^2$')
    loglog(plot_input_ary,error_ary_mat_vel(3,:),'rs','DisplayName','$\omega$, uniform')
    
    xlim([plot_input_ary(end) plot_input_ary(1)])
    
    ylabel('error'), xlabel('$\Delta t$')
    title(["Vel+TracIC: "+vel_type+"+"+IC_type,"IF-RK4+RK4SL; $Nx=$"+Nx])
    
    pplot(8,0.85,8)
    legend('Location','best','NumColumns',2)
    hold off
end

%%
figure(101)
pplot(8,0.78,8)
heatmap2d(IC_tracer_real,x_mesh,y_mesh); hold on
title("Initial Tracer Distrib $c(x,0)$")
xlabel("$x$"); ylabel("$y$")
pplot(8,0.78,8)

figure(102)
pplot(8,0.78,8)
heatmap2d(tracer_final,x_mesh,y_mesh); hold on
title("Final Tracer Distrib $c(x,T)$")
xlabel("$x$"); ylabel("$y$")
pplot(8,0.78,8)

%%
figure(100)
savefig("latex/figs/"+"NSAdv_conv_order_"+vel_type+"_"+IC_type)
figure(101)
savefig("latex/figs/"+"NSAdv_c_init_"+vel_type+"_"+IC_type)
figure(102)
savefig("latex/figs/"+"NSAdv_c_final_"+vel_type+"_"+IC_type)
