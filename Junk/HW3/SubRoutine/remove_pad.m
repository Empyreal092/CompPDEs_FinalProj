function [data_all_rm] = remove_pad(data_padall,Nx,Ny)

[My,Mx]=size(data_padall);
% tru_sz_x = (Mx+1)/3-1;
% tru_sz_y = (My+1)/3-1;
tru_sz_x = (Nx-1)/2;
tru_sz_y = (Ny-1)/2;

data_x_rm = [data_padall(:,1:tru_sz_x+1) data_padall(:,2*tru_sz_x+3:end)];
data_all_rm = [data_x_rm(1:tru_sz_y+1,:); data_x_rm(2*tru_sz_y+3:end,:)];

end

