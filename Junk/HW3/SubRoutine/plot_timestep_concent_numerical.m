omega_timestep_real = real(ifft2_n(concent_temp));

figure(1234)
hold on
omegaplot = pcolor(x_mesh,y_mesh,omega_timestep_real);
set(omegaplot, 'EdgeColor', 'none')
axis equal
colorbar
if T_timestep == dt
    pplot(12,0.8)
end