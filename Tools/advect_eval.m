function [jacob_four] = advect_eval(omega,concent)
arguments
    omega
    concent = omega;
end

global v0 L

[Ny,Nx] = size(omega);
k_ary = [0:(Nx-1)/2 -(Nx-1)/2:-1]*2*pi/L;
l_ary = [0:(Ny-1)/2 -(Ny-1)/2:-1]*2*pi/L;
[k_mesh,l_mesh] = meshgrid(k_ary,l_ary);

inv_laplace_mat = inv_laplace_make(k_ary,l_ary);

psi = inv_laplace_mat.*omega;
u_four = -1i*psi.*l_mesh;
v_four =  1i*psi.*k_mesh;
concent_px = 1i*concent.*k_mesh;
concent_py = 1i*concent.*l_mesh;

%%
% if if_anti_alias
    concent_px_pad_real = real(ifft2_n(zero_pad(concent_px)));
    concent_py_pad_real = real(ifft2_n(zero_pad(concent_py)));
    u_pad_real = v0+real(ifft2_n(zero_pad(u_four)));
    v_pad_real = v0+real(ifft2_n(zero_pad(v_four)));
    
    jacob_real = u_pad_real.*concent_px_pad_real+v_pad_real.*concent_py_pad_real;
    jacob_four = remove_pad(fft2_n(jacob_real),width(concent_px),height(concent_px));
% else
%     concent_px_pad_real = real(ifft2_n(concent_px));
%     concent_py_pad_real = real(ifft2_n(concent_py));
%     u_pad_real = v0+real(ifft2_n(u_four));
%     v_pad_real = v0+real(ifft2_n(v_four));
%     
%     jacob_real = u_pad_real.*concent_px_pad_real+v_pad_real.*concent_py_pad_real;
%     jacob_four = fft2_n(jacob_real);
% end
end

