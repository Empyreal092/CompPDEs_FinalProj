function [tracer_temp] = RK4SL_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,uh,vh,up,vp)
global L dt ext_sz finufft_interp

%%
if finufft_interp
    in_method_interp = "finufft";
else
    in_method_interp = "spline";
end

x1 = mod(x_mesh,L); y1 = mod(y_mesh,L);
% [u1,v1,~,~] = uvpomega_taylor(x1,y1,T_timestep+dt);
u1 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,up,x1,y1,in_method_interp);
v1 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vp,x1,y1,in_method_interp);

x2 = mod(x_mesh-dt/2*u1,L); y2 = mod(y_mesh-dt/2*v1,L);
% [u2,v2,~,~] = uvpomega_taylor(x2,y2,T_timestep+dt/2);
u2 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x2,y2,in_method_interp);
v2 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x2,y2,in_method_interp);

x3 = mod(x_mesh-dt/2*u2,L); y3 = mod(y_mesh-dt/2*v2,L);
% [u3,v3,~,~] = uvpomega_taylor(x3,y3,T_timestep+dt/2);
u3 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,uh,x3,y3,in_method_interp);
v3 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vh,x3,y3,in_method_interp);

x4 = mod(x_mesh-dt*u3,L); y4 = mod(y_mesh-dt*v3,L);
% [u4,v4,~,~] = uvpomega_taylor(x4,y4,T_timestep);
u4 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x4,y4,in_method_interp);
v4 = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x4,y4,in_method_interp);

um = (u1 + 2*u2 + 2*u3 + u4)/6;  vm = (v1 + 2*v2 + 2*v3 + v4)/6;

x_depart = mod((x_mesh - um*dt),L);
y_depart = mod((y_mesh - vm*dt),L);

tracer_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_depart,y_depart,"finufft");

end

