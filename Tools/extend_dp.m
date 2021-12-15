function [extended_data] = extend_dp(data,size)

left_idx = 1:size;
right_idx = width(data)-size+1:width(data);

data_left = data(:,left_idx);
data_right = data(:,right_idx);
data_top = data(left_idx,:);
data_bottom = data(right_idx,:);

data_tl = data(left_idx,left_idx);
data_tr = data(left_idx,right_idx);
data_bl = data(right_idx,left_idx);
data_br = data(right_idx,right_idx);

extended_data = [data_br data_bottom data_bl;
                 data_right data data_left;
                 data_tr data_top data_tl];

end

