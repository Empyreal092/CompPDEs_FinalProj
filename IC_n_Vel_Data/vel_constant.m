function [u,v,omega] = vel_constant(x,y,t)
global v0 Nx Ny

u = v0*ones(Nx,Ny);
v = v0*ones(Nx,Ny);
omega = zeros(Nx,Ny);

end

