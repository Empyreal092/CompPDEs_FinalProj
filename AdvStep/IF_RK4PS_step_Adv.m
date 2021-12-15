function [tracer_temp] = IF_RK4PS_step_Adv(tracer_temp,omega_n,omega_h,omega_p)
% global L dt ext_sz finufft_interp
global dt

omega_four_n = fft2_n(omega_n);
omega_four_h = fft2_n(omega_h);
omega_four_p = fft2_n(omega_p);

tracer_temp_four = fft2_n(tracer_temp);
%%
y1 = tracer_temp_four;
RHS_1 = -advect_eval(omega_four_n,y1);

y2 = tracer_temp_four+dt/2*RHS_1;
RHS_2 = -advect_eval(omega_four_h,y2);

y3 = tracer_temp_four+dt/2*RHS_2;
RHS_3 = -advect_eval(omega_four_h,y3);

y4 = tracer_temp_four+dt*RHS_3;
RHS_4 = -advect_eval(omega_four_p,y4);

RHSm = (RHS_1 + 2*RHS_2 + 2*RHS_3 + RHS_4)/6;
tracer_temp_four = tracer_temp_four + dt*RHSm;

tracer_temp = ifft2_n(tracer_temp_four);
end

