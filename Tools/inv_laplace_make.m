function inv_laplace_mat = inv_laplace_make(kx,ly)

[kkx,lly] = meshgrid(kx,ly);

inv_laplace_mat_neg = -kkx.^2-lly.^2;
inv_laplace_mat = inv_laplace_mat_neg.^-1;
inv_laplace_mat(1,1) = 0;

end

