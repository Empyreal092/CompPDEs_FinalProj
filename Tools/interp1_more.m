function [data_interp] = interp1_more(ext_sz,x_mesh_extend,data,x_depart,interp_method)
global L

Nx = length(x_depart);
if interp_method == "finufft"
    data_ft = fftshift(fft(data)/Nx);
    x_depart_ary = x_depart*2*pi/L;
    data_interp_ary = real(finufft1d2(x_depart_ary,1,1e-15,data_ft));
    data_interp = data_interp_ary';
else
    data_extend = extend_1p(data,ext_sz);
    data_interp = interp1(x_mesh_extend,data_extend,x_depart,interp_method);
end

end

