function [IC_concent_real] = IC_sinp1D(x_mesh,p)
global L

IC_concent_real = ( sin(pi*x_mesh/L) ).^p;

end

