function [IC_omega_real] = IC_3vort(x_mesh,y_mesh)
global L
% 
vor_1 = @(x_mesh,y_mesh) exp( -5*((x_mesh-pi).^2+(y_mesh-3*pi/4).^2) );
vor_2 = @(x_mesh,y_mesh) exp( -5*((x_mesh-pi).^2+(y_mesh-5*pi/4).^2) );
vor_3 = @(x_mesh,y_mesh) -1/2*exp( -5/2*((x_mesh-5*pi/4).^2+(y_mesh-5*pi/4).^2) );
% vor_1 = @(x_mesh,y_mesh) exp( -1*((x_mesh-pi).^2+(y_mesh-3*pi/4).^2) );
% vor_2 = @(x_mesh,y_mesh) exp( -1*((x_mesh-pi).^2+(y_mesh-5*pi/4).^2) );
% vor_3 = @(x_mesh,y_mesh) -1/2*exp( -1/2*((x_mesh-5*pi/4).^2+(y_mesh-5*pi/4).^2) );

IC_omega = zeros(size(x_mesh));
% for i = -0:0
%     for j = -0:0
for i = -10:10
    for j = -10:10
        x_mesh_move = (x_mesh+i*L);
        y_mesh_move = (y_mesh+j*L);
        IC_omega = IC_omega + vor_1(x_mesh_move,y_mesh_move)...
            +vor_2(x_mesh_move,y_mesh_move)+vor_3(x_mesh_move,y_mesh_move);
    end
end

IC_omega_real = IC_omega;

%%
% omegaplot = pcolor(x_mesh,y_mesh,IC_omega);
% set(omegaplot, 'EdgeColor', 'none')
% axis equal
% colorbar

end

