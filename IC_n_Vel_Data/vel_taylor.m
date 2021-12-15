function [u,v,omega] = vel_taylor(x,y,t)
global L v0

% v0 = 1;
mu = 0;

exp_fac = exp( -8*pi^2*mu*t/L^2 );
x_val = 2*pi*(x-v0*t)/L;
y_val = 2*pi*(y-v0*t)/L;

u = v0 - 2*exp_fac*cos(x_val).*sin(y_val);
v = v0 + 2*exp_fac*sin(x_val).*cos(y_val);
% p = -exp_fac^2*( cos(2*x_val)+cos(2*y_val) );
% 
omega = 2*exp_fac*4*pi*cos(x_val).*cos(y_val)/L;

end

