function [tracer_temp] = Trap_step_Adv(tracer_temp,x_mesh,y_mesh,x_mesh_extend,y_mesh_extend,un,vn,um,vm)
global L dt ext_sz finufft_interp

iter_count = 0;

up = 2*un-um; vp = 2*vn-vm;
un_iter = un; vn_iter = vn;
while iter_count < 1 % only one iteration
    x_depart_iter = mod(x_mesh - (un_iter+up)/2*dt,L);
    y_depart_iter = mod(y_mesh - (vn_iter+vp)/2*dt,L);
    if finufft_interp
        un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_depart_iter,y_depart_iter,"finufft");
        vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_depart_iter,y_depart_iter,"finufft");
    else
        un_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,un,x_depart_iter,y_depart_iter,"linear");
        vn_iter = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,vn,x_depart_iter,y_depart_iter,"linear");
    end
    iter_count = iter_count+1;
end
x_depart = mod((x_mesh - (un_iter+up)/2*dt),L);
y_depart = mod((y_mesh - (vn_iter+vp)/2*dt),L);

if finufft_interp
    tracer_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_depart,y_depart,"finufft");
else
    tracer_temp = interp2_more(ext_sz,x_mesh_extend,y_mesh_extend,tracer_temp,x_depart,y_depart,"cubic");
end

end

