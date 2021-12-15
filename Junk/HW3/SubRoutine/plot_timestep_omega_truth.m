[~,~,~,omega_truth] = taylor_exact(x_mesh,y_mesh,T_timestep);

figure(1235)
hold on
omegaplot = pcolor(x_mesh,y_mesh,omega_truth);
set(omegaplot, 'EdgeColor', 'none')
axis equal
colorbar
if T_timestep == dt
    pplot(12,0.8)
end