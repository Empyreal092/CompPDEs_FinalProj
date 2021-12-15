function [IC_omega_four] = IC_omega_Taylor(x_mesh,y_mesh)

[~,~,~,IC_omega] = taylor_exact(x_mesh,y_mesh,0);
IC_omega_four = fft2_n(IC_omega);

end

