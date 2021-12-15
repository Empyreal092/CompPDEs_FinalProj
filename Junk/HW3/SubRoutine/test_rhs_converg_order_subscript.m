Mx = 2^6+1; My = Mx;
x_dense = 0:L/Mx:L-L/Mx; y_dense = 0:L/My:L-L/My;
[x_mesh_d,y_mesh_d] = meshgrid(x_dense,y_dense);
k_dense = [0:(Mx-1)/2 -(Mx-1)/2:-1]*2*pi/L; l_dense = [0:(My-1)/2 -(My-1)/2:-1]*2*pi/L;
[k_mesh_d,l_mesh_d] = meshgrid(k_dense,l_dense);
k2l2_d = -k_mesh_d.^2-l_mesh_d.^2;
A_d = mu*k2l2_d;

IC_omega_dense = IC_omega_Vortices(x_mesh_d,y_mesh_d);
adv_truth = advect_eval(IC_omega_dense);
diff_truth = A_d.*IC_omega_dense;

rhs_Nxy_ary = [1:2:Mx-2];
error_rhs_adv_ary_l1 = []; error_rhs_adv_ary_l2 = []; error_rhs_adv_ary_inf = [];
error_rhs_diff_ary_l1 = []; error_rhs_diff_ary_l2 = []; error_rhs_diff_ary_inf = [];
for rhs_Nx = rhs_Nxy_ary
    rhs_Ny = rhs_Nx;
    [rhs_x_mesh,rhs_y_mesh] = meshgrid(0:L/rhs_Nx:L-L/rhs_Nx,0:L/rhs_Ny:L-L/rhs_Ny);
    k_rhs = [0:(rhs_Nx-1)/2 -(rhs_Nx-1)/2:-1]*2*pi/L; l_rhs = [0:(rhs_Ny-1)/2 -(rhs_Ny-1)/2:-1]*2*pi/L;
    [k_mesh_rhs,l_mesh_rhs] = meshgrid(k_rhs,l_rhs);
    k2l2_rhs = -k_mesh_rhs.^2-l_mesh_rhs.^2;
    A_rhs = mu*k2l2_rhs;
    
    IC_omega_test = IC_omega_Vortices(rhs_x_mesh,rhs_y_mesh);
    adv_test = advect_eval(IC_omega_test);
    adv_test_interpn = zero_pad(adv_test,Mx,My);
    diff_test = A_rhs.*IC_omega_test;
    diff_test_interpn = zero_pad(diff_test,Mx,My);
    
    error_rhs_adv_l1 = grid_func_norm(ifft2_n(adv_test_interpn),ifft2_n(adv_truth),"1");
    error_rhs_adv_l2 = grid_func_norm(ifft2_n(adv_test_interpn),ifft2_n(adv_truth),"2");
    error_rhs_adv_inf = grid_func_norm(ifft2_n(adv_test_interpn),ifft2_n(adv_truth),"inf");
    error_rhs_adv_ary_l1 = [error_rhs_adv_ary_l1 error_rhs_adv_l1];
    error_rhs_adv_ary_l2 = [error_rhs_adv_ary_l2 error_rhs_adv_l2];
    error_rhs_adv_ary_inf = [error_rhs_adv_ary_inf error_rhs_adv_inf];
    error_rhs_diff_l1 = grid_func_norm(ifft2_n(diff_test_interpn),ifft2_n(diff_truth),"1");
    error_rhs_diff_l2 = grid_func_norm(ifft2_n(diff_test_interpn),ifft2_n(diff_truth),"2");
    error_rhs_diff_inf = grid_func_norm(ifft2_n(diff_test_interpn),ifft2_n(diff_truth),"inf");
    error_rhs_diff_ary_l1 = [error_rhs_diff_ary_l1 error_rhs_diff_l1];
    error_rhs_diff_ary_l2 = [error_rhs_diff_ary_l2 error_rhs_diff_l2];
    error_rhs_diff_ary_inf = [error_rhs_diff_ary_inf error_rhs_diff_inf];
end

if if_anti_alias
    figure(99)
else
    figure(98)
end
semilogy(rhs_Nxy_ary,error_rhs_adv_ary_l1,'ro','DisplayName','$J(\psi,\omega), \ell^1$'); hold on
semilogy(rhs_Nxy_ary,error_rhs_adv_ary_l2,'r^','DisplayName','$J(\psi,\omega), \ell^2$')
semilogy(rhs_Nxy_ary,error_rhs_adv_ary_inf,'rs','DisplayName','$J(\psi,\omega),$ uniform')
semilogy(rhs_Nxy_ary,error_rhs_diff_ary_l1,'bo','DisplayName','$\mu\nabla^2\omega, \ell^1$')
semilogy(rhs_Nxy_ary,error_rhs_diff_ary_l2,'b^','DisplayName','$\mu\nabla^2\omega, \ell^2$')
semilogy(rhs_Nxy_ary,error_rhs_diff_ary_inf,'bs','DisplayName','$\mu\nabla^2\omega,$ uniform')
legend('Location','best','NumColumns',2)
ylabel('error'), xlabel('$Nx,Ny$')
if if_anti_alias
    title(["Errors in Evaluating RHSs for the 3 Vortices IC","Anti-Aliasing: On"])
else
    title(["Errors in Evaluating RHSs for the 3 Vortices IC","Anti-Aliasing: Off"])
end
ylimm = ylim;
ylim([1e-20 2e0])
pplot(8,0.9,8)
hold off
