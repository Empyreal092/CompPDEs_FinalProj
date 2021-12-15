function [data_padall] = zero_pad(data,Mx,My)

arguments
    data
    Mx = (width(data)+1)/2*3-1;
    My = (height(data)+1)/2*3-1;
end

[Ny,Nx]=size(data);
cut_x = (Nx+1)/2;
cut_y = (Ny+1)/2;

% padsz_x = cut_x;
% padsz_y = cut_y;
padsz_x = Mx-Nx;
padsz_y = My-Ny;

data_padx = [data(:,1:cut_x) zeros(Ny,padsz_x) data(:,cut_x+1:end)];
data_padall = [data_padx(1:cut_y,:);zeros(padsz_y,Nx+padsz_x);data_padx(cut_y+1:end,:)];

end

