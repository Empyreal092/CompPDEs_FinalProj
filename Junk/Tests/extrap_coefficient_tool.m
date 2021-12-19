x = -1/2;

x0 = 0; x1 = -1; x2 = -2; x3 = -3;

l0 = @(x) (x-x1)/(x0-x1)*(x-x2)/(x0-x2)*(x-x3)/(x0-x3);
l1 = @(x) (x-x0)/(x1-x0)*(x-x2)/(x1-x2)*(x-x3)/(x1-x3);
l2 = @(x) (x-x0)/(x2-x0)*(x-x1)/(x2-x1)*(x-x3)/(x2-x3);
l3 = @(x) (x-x0)/(x3-x0)*(x-x1)/(x3-x1)*(x-x2)/(x3-x2);

% a0 = l0(x)*16;
% a1 = l1(x)*16;
% a2 = l2(x)*16;
% a3 = l3(x)*16;
const = 16;

a0 = l0(x)*const;
a1 = l1(x)*const;
a2 = l2(x)*const;
a3 = l3(x)*const;

disp(a0+" "+a1+" "+a2+" "+a3)

%%
% g = @(x) x^3+1-(x+10)^2;
% 
% g0 = g(0); g1 = g(-1); g2 = g(-2); g3 = g(-3);
% 
% gh = 1/16*(35*g0-35*g1+21*g2-5*g3)
% g(1/2)
% 
% gp = (4*g0-6*g1+4*g2-g3)
% g(1)
% %         omega_temp_h = 1/16*(35*omega_temp-35*omega_temp_prev+21*omega_temp_prev2-5*omega_temp_prev3);
% %         omega_temp_p = (4*omega_temp-6*omega_temp_prev+4*omega_temp_prev2-omega_temp_prev3);