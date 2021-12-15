function [omega_temp] = IF_RK4PS_step_NS(omega_temp)
global dt

omega_four_n = fft2_n(omega_temp);
% omega_four_h = fft2_n(omega_h);
% omega_four_p = fft2_n(omega_p);

% tracer_temp_four = fft2_n(tracer_temp);
%%
y1 = omega_four_n;
RHS_1 = -advect_eval(y1);

y2 = omega_four_n+dt/2*RHS_1;
RHS_2 = -advect_eval(y2);

y3 = omega_four_n+dt/2*RHS_2;
RHS_3 = -advect_eval(y3);

y4 = omega_four_n+dt*RHS_3;
RHS_4 = -advect_eval(y4);

RHSm = (RHS_1 + 2*RHS_2 + 2*RHS_3 + RHS_4)/6;
omega_four_n = omega_four_n + dt*RHSm;

omega_temp = ifft2_n(omega_four_n);

end

