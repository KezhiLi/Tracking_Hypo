function texture = ske2tex(xy, width, Y_k)
% 
% 
% 
% 
% 
% 
% 
% 



num_pt = size(xy,1);

Y_k_gray = 255 - rgb2gray(Y_k);

% neglect first 2 and last 2 points
neglect_pt = 2; 

num_step = (num_pt-2*neglect_pt)*2-1;
pt_on_ske = zeros(num_step,2);

for ii = 1:num_pt-neglect_pt*2;
    pt_on_ske(ii*2-1,:) = xy(ii+neglect_pt,:);
end
for ii = 2:2:size(pt_on_ske,1);
    pt_on_ske(ii,:) = (pt_on_ske(ii-1,:)+pt_on_ske(ii+1,:))/2;
end

stat_pt = zeros(num_step,2);


for jj = 1:num_step;
    mask_cir = create_cir(pt_on_ske(jj,2), pt_on_ske(jj,1), width, size(Y_k_gray));
    cir_full = double(Y_k_gray).* mask_cir;
    cir = cir_full(cir_full>0);
    stat_pt(jj,1) = mean2(cir);
    stat_pt(jj,2) = std2(cir);
end

texture = stat_pt;



