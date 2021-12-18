function [IC_concent_real] = IC_step(x_mesh)

IC_concent_real = zeros(1,length(x_mesh));
for i = 1:length(x_mesh)
    if x_mesh(i) > 1/4 && x_mesh(i) < 1/2
        IC_concent_real(i) = 1;
    end
end

end

