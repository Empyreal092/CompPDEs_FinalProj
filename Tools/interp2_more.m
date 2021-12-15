function [data_interp] = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,data,x_depart,y_depart,interp_method)
global L

[Nx,Ny] = size(x_depart);
if interp_method == "finufft"
    data_ft = fftshift(fft2_n(data));
    x_depart_ary = reshape(x_depart,1,numel(x_depart))*2*pi/L;
    y_depart_ary = reshape(y_depart,1,numel(y_depart))*2*pi/L;
    data_interp_ary = real(finufft2d2(y_depart_ary,x_depart_ary,1,1e-15,data_ft));
    data_interp = reshape(data_interp_ary,Nx,Ny);
else
    data_extend = extend_dp(data,ext_sz);
    data_interp = interp2(x_mesh_extend,y_mesh_extend,data_extend,x_depart,y_depart,interp_method);
end

end

