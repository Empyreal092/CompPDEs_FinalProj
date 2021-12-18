function [extended_data] = extend_1p(data,size)

left_idx = 1:size;
right_idx = width(data)-size+1:width(data);

data_left = data(left_idx);
data_right = data(right_idx);

extended_data = [data_right data data_left];

end

