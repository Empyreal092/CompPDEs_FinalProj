close all
clear all

error_ary = [];

N_ary = 2.^[1:10];
for N = N_ary

x = 0:1/N:1;
v = @(x) sin(10*x)+x;

xq = 1+1/N;

vq = interp1(x,v(x),xq,'spline','extrap');

error_ary = [error_ary abs(vq-v(xq))];

end

%%
figure
loglog(N_ary,N_ary.^-4,'k'); hold on
loglog(N_ary,error_ary,'o')

%%
% figure
% plot([x xq],v([x xq])); hold on
% plot([x(end) xq],[v(end) vq])
