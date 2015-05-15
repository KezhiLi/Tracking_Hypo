function blk = cropBody(xy, tail_crop, num_crop, width, Y_k)
% 
% 
%
% 
% 
%  num_crop < len_seg - 1, tail_crop + num_crop < len_seg

Npix_h = size(Y_k, 1);
Npix_w = size(Y_k, 2);

blank = logical(zeros(Npix_h,0)*zeros(0,Npix_w));

%len_xy = size(xy, 1);
%xy_consid = xy(tail_crop:tail_crop+num_crop,:);

xy_con_plus = xy(tail_crop:tail_crop+num_crop,:) + xy(tail_crop+1:tail_crop+num_crop+1,:);

centr_pt = xy_con_plus/2;

for ii = 1: num_crop;
    blk_ii = blank;
    pt_ii = centr_pt(ii,:);
    blk_ii(round(pt_ii(1)-width):round(pt_ii(1)+width),round(pt_ii(2)-width):round(pt_ii(2)+width)) = 1;
    blk{ii} = Y_k .* blk_ii;
end








