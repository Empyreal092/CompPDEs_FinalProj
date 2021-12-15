function [IC_concent_four] = IC_concent(x_mesh,y_mesh,p)
global L

IC_concent_real = ( sin(pi*x_mesh/L).*sin(pi*y_mesh/L) ).^p;
IC_concent_four = fft2_n(IC_concent_real);

end

