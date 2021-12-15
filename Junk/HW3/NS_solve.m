clear all
% close all
addpath(genpath('SubRoutine'))  

global L v0 mu if_anti_alias

%%
IC_type = "Taylor";
% IC_type = "Vortices";

% time_step_method = "SBDF2";
time_step_method = "IF-RK4";

if_anti_alias = true;

plot_timestep_numr = false;
plot_timestep_truth = false;

test_time_converg_order = false;
test_rhs_converg_order = true;

%%
switch IC_type
    case "Taylor"
        L = 1;
        v0 = 1;
    case "Vortices"
        L = 2*pi;
        v0 = 1;
end
T = 0.25*1e0;
mu = 0.05;
% mu = 0;

error_ary_l1_uv = []; error_ary_l2_uv = []; error_ary_inf_uv = [];
error_ary_l1_omega = []; error_ary_l2_omega = []; error_ary_inf_omega = [];

%%
if test_time_converg_order
    Nt_pow = [1:7];
else
    Nt_pow = 3;
end
Nt_ary = 2.^Nt_pow; plot_input_ary = Nt_ary;

if test_time_converg_order
    Mx = 2^10+1; My = Mx;
    omega_final_save = -1e15*ones(Mx,My);
end

for Nt = Nt_ary
    Nx = 2^5+1;
    Ny = Nx;
    %     Nt = 2^10;
    dt = T/Nt;
    
    x_ary = 0:L/Nx:L-L/Nx;
    y_ary = 0:L/Ny:L-L/Ny;
    [x_mesh,y_mesh] = meshgrid(x_ary,y_ary);
    k_ary = [0:(Nx-1)/2 -(Nx-1)/2:-1]*2*pi/L;
    l_ary = [0:(Ny-1)/2 -(Ny-1)/2:-1]*2*pi/L;
    [k_mesh,l_mesh] = meshgrid(k_ary,l_ary);
    
    switch IC_type
        case "Taylor"
            IC_omega_four = IC_omega_Taylor(x_mesh,y_mesh);
        case "Vortices"
            IC_omega_four = IC_omega_Vortices(x_mesh,y_mesh);
    end
    
    %%
    if test_rhs_converg_order && IC_type == "Vortices"
        test_rhs_converg_order_subscript
    end
    
    %%
    T_timestep = 0;
    k2l2 = -k_mesh.^2-l_mesh.^2;
    A = mu*k2l2;
    
    switch time_step_method
        case "SBDF2"
            omega_temp = IC_omega_four;
            omega_2prevsteps = -1e15*ones(Nx,Nx,2); % place to store the last 2 solutions
            omega_2prevsteps(:,:,2) = IC_omega_four; % the 2nd column is the most recent sol
            B_2prevsteps = -1e15*ones(Nx,Nx,2); % place to store the last 2 solutions
            B_2prevsteps(:,:,2) = -advect_eval(IC_omega_four); % the 2nd column is the most recent sol

            while T_timestep < T-dt/2
                B = B_2prevsteps(:,:,2);
                if omega_2prevsteps(1,1) < -1e9 % if it is the first step, use Euler for advection, and Crank-Nicolson for diffusion instead
                    omega_temp = ( omega_temp + dt/2*A.*omega_temp + dt*B )./(ones(Nx,Nx)-dt/2*A);
                else % otherwise, use SBDF2
                    omega_temp = ( 4/3*omega_2prevsteps(:,:,2)-1/3*omega_2prevsteps(:,:,1) + 2/3*dt*(2*B_2prevsteps(:,:,2)-B_2prevsteps(:,:,1)) )...
                        ./(ones(Nx,Nx)-2*dt/3*A);
                end
                
                % update the history array
                omega_2prevsteps(:,:,1) = omega_2prevsteps(:,:,2);
                omega_2prevsteps(:,:,2) = omega_temp;
                B_2prevsteps(:,:,1) = B_2prevsteps(:,:,2);
                B_2prevsteps(:,:,2) = -advect_eval(omega_temp);
                
                T_timestep = T_timestep+dt;
                %%
                if plot_timestep_numr
                    plot_timestep_omega_numerical
                end
                %%
                if plot_timestep_truth
                    plot_timestep_omega_truth
                end
            end
            omega_final = omega_temp;
        case "IF-RK4"
            omega_temp = IC_omega_four;
            
            dth = dt/2;
            itfmh = exp(A*dth); itfm1 = exp(A*dt); itfm0 = exp(A*0);
            
            while T_timestep < T-dt/2
                y1 = omega_temp;
                RHS_1 = -advect_eval(y1);
                y2 = itfmh.*omega_temp+dth*itfmh.*RHS_1;
                RHS_2 = -advect_eval(y2);
                y3 = itfmh.*omega_temp+dth*itfm0.*RHS_2;
                RHS_3 = -advect_eval(y3);
                y4 = itfm1.*omega_temp+dt *itfmh.*RHS_3;
                RHS_4 = -advect_eval(y4);
                
                RHSm = (itfm1.*RHS_1 + 2*itfmh.*RHS_2 + 2*itfmh.*RHS_3 + itfm0.*RHS_4)/6;
                omega_temp = itfm1.*omega_temp + dt*RHSm;
                
                T_timestep = T_timestep+dt;
                %%
                if plot_timestep_numr
                    plot_timestep_omega_numerical
                end
                %%
                if plot_timestep_truth
                    plot_timestep_omega_truth
                end
            end
            omega_final = omega_temp;
    end
    %%
    if test_time_converg_order && IC_type == "Taylor"
        test_time_converg_order_taylor_subscipt
    end
    if test_time_converg_order && IC_type == "Vortices"
        test_time_converg_order_vortices_subscipt
        omega_final_save = zero_pad(omega_final,Mx,My);
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
                loglog(plot_input_ary,50*plot_input_ary.^-2,'k','DisplayName','2nd Order'); hold on
            else
                loglog(plot_input_ary,1/5*plot_input_ary.^-2,'k','DisplayName','2nd Order'); hold on
            end
        case "IF-RK4"
            if IC_type == "Taylor"
                loglog(plot_input_ary,5*plot_input_ary.^-4,'k','DisplayName','4th Order'); hold on
            else
                loglog(plot_input_ary,1/5e1*plot_input_ary.^-4,'k','DisplayName','4th Order'); hold on
            end
    end
    
    loglog(plot_input_ary,error_ary_l1_uv,'ro','DisplayName','$(u,v),\ell^1$')
    loglog(plot_input_ary,error_ary_l2_uv,'r^','DisplayName','$(u,v),\ell^2$')
    loglog(plot_input_ary,error_ary_inf_uv,'rs','DisplayName','$(u,v)$, uniform')
    
    loglog(plot_input_ary,error_ary_l1_omega,'bo','DisplayName','$\omega,\ell^1$')
    loglog(plot_input_ary,error_ary_l2_omega,'b^','DisplayName','$\omega,\ell^2$')
    loglog(plot_input_ary,error_ary_inf_omega,'bs','DisplayName','$\omega$, uniform')
    
    ylabel('error'), xlabel('$Nt$')
    title("IC: "+IC_type+"; "+time_step_method)
    xlim([plot_input_ary(1) plot_input_ary(end)])
    ylimm = ylim;
    ylim([ylimm(1)/100 ylimm(end)])
    legend('Location','best','NumColumns',2)
    pplot(8,1,8)
end

