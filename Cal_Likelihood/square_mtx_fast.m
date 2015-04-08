function dis_mtx = square_mtx_fast(dis_mtx, pt, len)
%
% 
% 
% 
% 

Y_1 = size(dis_mtx, 1);
Y_2 = size(dis_mtx, 2);


left = round(max(pt(1)-len,1));
right =round(min(pt(1)+len,Y_2));
bottom = round(max(pt(2)-len,1));
top = round(min(pt(2)+len,Y_1));

dis_mtx(bottom:top,left:right) = 2;


