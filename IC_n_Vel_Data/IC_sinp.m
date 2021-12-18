function [IC_concent_real] = IC_sinp(x_mesh,y_mesh,p)
global L

IC_concent_real = ( sin(pi*x_mesh/L).*sin(pi*y_mesh/L) ).^p;

end

