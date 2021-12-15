close all

L = 1;
N_ary = round(1.5.^[4:15]);
error_ary = [];

for N = N_ary
% N = 100;
h = L/N;

x_ary = 0:h:L-h;
y_ary = 0:h:L-h;
[x_mesh,y_mesh] = meshgrid(x_ary,y_ary);

% func = @(x,y) sin(2*pi*x).*sin(2*pi*y);
func = @(x,y) sin(2*pi*x);

func_val = func(x_mesh,y_mesh);

%% 
% figure(1234)
% hold on
% omegaplot = pcolor(x_mesh,y_mesh,func_val);
% set(omegaplot, 'EdgeColor', 'none')
% axis equal
% colorbar
% pplot(12,0.8)

%% 
x_q = sqrt(2)/2; y_q = 0;
interp_val = interp2(x_mesh,y_mesh,func_val,x_q,y_q,'cubic');
exact_val = func(x_q,y_q);

error_ary = [error_ary abs(interp_val-exact_val)];
end

loglog(N_ary,error_ary,'o'); hold on
loglog(N_ary,N_ary.^-2)
loglog(N_ary,N_ary.^-3)
loglog(N_ary,N_ary.^-4)
pplot(12,0.8)

