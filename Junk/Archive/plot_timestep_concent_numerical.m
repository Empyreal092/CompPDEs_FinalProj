figure(1234)
hold on
omegaplot = pcolor();
set(omegaplot, 'EdgeColor', 'none')
axis equal
colorbar
if T_timestep == dt
    pplot(12,0.8)
end