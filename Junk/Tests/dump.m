%%
%     uh = 1/8*(15*un-10*um+3*um2); vh = 1/8*(15*vn-10*vm+3*vm2);
%     up = (3*un-3*um+um2);     vp = (3*vn-3*vm+vm2);
% 


%     uh = zeros(Nx,Ny); vh = zeros(Nx,Ny); up = zeros(Nx,Ny); vp = zeros(Nx,Ny); 
%     for i = 1:Ny
%         for j = 1:Nx
%             t_ary = [0 -1 -2 -3];
%             u_ary = [un(i,j) um(i,j) um2(i,j) um3(i,j)];
%             v_ary = [vn(i,j) vm(i,j) vm2(i,j) vm3(i,j)];
%             
%             tq = [0.5 1];
%             vqu = interp1(t_ary,u_ary,tq,'spline','extrap'); uh(i,j) = vqu(1); up(i,j) = vqu(2);
%             vqv = interp1(t_ary,v_ary,tq,'spline','extrap'); vh(i,j) = vqv(1); vp(i,j) = vqv(2);
%         end
%     end


    %% Extrapolate Omega
%     omega_temp_h = 1/16*(35*omega_temp-35*omega_temp_prev+21*omega_temp_prev2-5*omega_temp_prev3);
%     omega_temp_p = (4*omega_temp-6*omega_temp_prev+4*omega_temp_prev2-omega_temp_prev3);
%     
%     omega_temp_h = 1/8*(15*omega_temp-10*omega_temp_prev+3*omega_temp_prev2);
%     omega_temp_p = (3*omega_temp-3*omega_temp_prev+1*omega_temp_prev2);
%     
%     [~,~,~,omega_temp_h1] = uvpomega_taylor(x_mesh,y_mesh,T_timestep+dt/2);
%     [~,~,~,omega_temp_p1] = uvpomega_taylor(x_mesh,y_mesh,T_timestep+dt);
%     
%     [un,vn] = omega_2_uv(omega_temp);
%     [uh,vh] = omega_2_uv(omega_temp_h);
%     [up,vp] = omega_2_uv(omega_temp_p);
    