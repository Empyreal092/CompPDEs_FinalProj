function [tracer_temp] = Euler_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn)
global L dt ext_sz finufft_interp

x_depart = mod(x_mesh - un*dt , L);
y_depart = mod(y_mesh - vn*dt , L);

if finufft_interp
    tracer_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_depart,y_depart,"finufft");
else
    tracer_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_depart,y_depart,"linear");
end

end

