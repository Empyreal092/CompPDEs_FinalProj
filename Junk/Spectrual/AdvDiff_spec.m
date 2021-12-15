clear all
close all
addpath(genpath('SubRoutine'))  

global L v0 mu if_anti_alias

%%
IC_type = "Taylor";
% IC_type = "Constant"; 

% time_step_method = "SBDF2";
time_step_method = "IF-RK4";

if_anti_alias = true;

plot_timestep_numr = true;
plot_timestep_truth = false;

test_time_converg_order = false;
test_rhs_converg_order = false;

%%
L = 1;
v0 = 1;
switch IC_type
    case "Taylor"
        mu = 0;
        T = 0.1;
    case "Constant"
        mu = 0;
        T = 1;
end

error_ary_l1_uv = []; error_ary_l2_uv = []; error_ary_inf_uv = [];

%%
if test_time_converg_order
    Nt_pow = [5:0.2:7];
else
    Nt_pow = 6;
end
Nt_ary = round(2.^Nt_pow); plot_input_ary = Nt_ary;

if test_time_converg_order
    Mx = 2^10+1; My = Mx;
    omega_final_save = -1e15*ones(Mx,My);
end
for Nt = Nt_ary
    Nx = 2^6+1;
    Ny = Nx;
    %     Nt = 2^10;
    dt = T/Nt;
    
    x_ary = 0:L/Nx:L-L/Nx;
    y_ary = 0:L/Ny:L-L/Ny;
    [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
    k_ary = [0:(Nx-1)/2 -(Nx-1)/2:-1]*2*pi/L;
    l_ary = [0:(Ny-1)/2 -(Ny-1)/2:-1]*2*pi/L;
    [k_mesh,l_mesh] = meshgrid(k_ary,l_ary);
    
%     switch IC_type
%         case "Taylor"
%             [u,~,~,~] = taylor_exact(x_mesh,y_mesh,0);
%             IC_concen_four = fft2_n(u);
%         case "Constant"
%             p = 20;
%             IC_concen_four = IC_concent(x_mesh,y_mesh,p);
%     end
    p = 2;
    IC_concen_four = IC_concent(x_mesh,y_mesh,p);
    
    %%
    T_timestep = 0;
    k2l2 = -k_mesh.^2-l_mesh.^2;
    A = mu*k2l2;
    
    switch time_step_method
        case "IF-RK4"
            concent_temp = IC_concen_four;
            
            dth = dt/2;
            itfmh = exp(A*dth); itfm1 = exp(A*dt); itfm0 = exp(A*0);
            
            while T_timestep < T-dt/2
                switch IC_type
                    case "Taylor"
                        [omega_data_four,~] = concent_Taylor_vel_data(x_mesh,y_mesh,T_timestep); forcing_four = 0;
                        [omega_data_four_h,~] = concent_Taylor_vel_data(x_mesh,y_mesh,T_timestep+dth); forcing_four_h = 0;
                        [omega_data_four_1,~] = concent_Taylor_vel_data(x_mesh,y_mesh,T_timestep+dt); forcing_four_1 = 0;
                end
                
                y1 = concent_temp;
                RHS_1 = -advect_eval(omega_data_four,y1)-forcing_four;
                y2 = itfmh.*concent_temp+dth*itfmh.*RHS_1;
                RHS_2 = -advect_eval(omega_data_four_h,y2)-forcing_four_h;
                y3 = itfmh.*concent_temp+dth*itfm0.*RHS_2;
                RHS_3 = -advect_eval(omega_data_four_h,y3)-forcing_four_h;
                y4 = itfm1.*concent_temp+dt *itfmh.*RHS_3;
                RHS_4 = -advect_eval(omega_data_four_1,y4)-forcing_four_1;
                
                RHSm = (itfm1.*RHS_1 + 2*itfmh.*RHS_2 + 2*itfmh.*RHS_3 + itfm0.*RHS_4)/6;
                concent_temp = itfm1.*concent_temp + dt*RHSm;
                
                T_timestep = T_timestep+dt;
                %%
                if plot_timestep_numr
                    plot_timestep_concent_numerical
                end
                %%
                if plot_timestep_truth
                    plot_timestep_concent_truth
                end
            end
            concent_final = concent_temp;
    end
    %%
    if test_time_converg_order && IC_type == "Taylor"
        test_time_converg_order_concent_taylor_subscipt
    end
    if test_time_converg_order && IC_type == "Constant"
        test_time_converg_order_concent_const_subscipt
    end
end

%%
if test_time_converg_order
    figure(100)
    
    switch IC_type
        case "Taylor"
            plot_input_ary = Nt_ary;
        case "Vortices"
            plot_input_ary = Nt_ary(2:end);
            error_ary_l1_uv = error_ary_l1_uv(2:end);
            error_ary_l2_uv = error_ary_l2_uv(2:end);
            error_ary_inf_uv = error_ary_inf_uv(2:end);
            error_ary_l1_omega = error_ary_l1_omega(2:end);
            error_ary_l2_omega = error_ary_l2_omega(2:end);
            error_ary_inf_omega = error_ary_inf_omega(2:end);
    end

    switch time_step_method
        case "SBDF2"
            if IC_type == "Taylor"
                loglog(plot_input_ary,5*plot_input_ary.^-2,'k','DisplayName','2nd Order'); hold on
            else
%                 loglog(plot_input_ary,500*plot_input_ary.^-2,'k','DisplayName','2nd Order'); hold on
                loglog(plot_input_ary,5e5*plot_input_ary.^-2,'k','DisplayName','2nd Order'); hold on
            end
        case "IF-RK4"
            if IC_type == "Taylor"
                loglog(plot_input_ary,5*plot_input_ary.^-4,'k','DisplayName','4th Order'); hold on
            else
                loglog(plot_input_ary,5e1*plot_input_ary.^-4,'k','DisplayName','4th Order'); hold on
%                 loglog(plot_input_ary,5e8*plot_input_ary.^-4,'k','DisplayName','4th Order'); hold on
            end
    end
    
    loglog(plot_input_ary,error_ary_l1_uv,'ro','DisplayName','$u,\ell^1$')
    loglog(plot_input_ary,error_ary_l2_uv,'r^','DisplayName','$u,\ell^2$')
    loglog(plot_input_ary,error_ary_inf_uv,'rs','DisplayName','$u$, uniform')
    
    ylabel('error'), xlabel('$Nt$')
    switch IC_type
        case "Taylor"
            title("IC: "+IC_type+"; "+time_step_method)
        case "Constant"
            title("IC: "+IC_type+"; "+time_step_method+"; $p=$"+p+"; $Nx=Ny=$"+Nx)
    end
    xlim([plot_input_ary(1) plot_input_ary(end)])
    ylimm = ylim;
%     ylim([ylimm(1) 1e1])
    legend('Location','best','NumColumns',1)
    pplot(8,0.8,8)
    hold off
end

