function ske_pred = re_interp(len_ske_pred, len_ske_old, ske_pred_1)
% 
% 
% 
% 


pp = 1:len_ske_pred;
qq = 1:((len_ske_pred-1)/(len_ske_old-1)):len_ske_pred;
ske_pred_x0 = interp1(pp,ske_pred_1(:,1),qq,'spline');
ske_pred_y0 = interp1(pp,ske_pred_1(:,2),qq,'spline');
ske_pred_2 = [ske_pred_x0',ske_pred_y0'];

qq2 = zeros(len_ske_old,1);
[length_ske_pred, dis_ske] = pt_len(ske_pred_2);
reci_dis_ske = [0;1./dis_ske];
sum_reci_dis_ske = sum(reci_dis_ske);

% in case len_ske_old > length(reci_dis_ske)
len_reci_dis_ske = length(reci_dis_ske);
for q = 1:len_ske_old;
    qq2(q) = (len_ske_old-1)*(sum(reci_dis_ske(1:min(q, len_reci_dis_ske))))/sum_reci_dis_ske+1;
end

ske_pred_x1 = interp1(1:len_ske_old,ske_pred_x0,qq2,'spline');
ske_pred_y1 = interp1(1:len_ske_old,ske_pred_y0,qq2,'spline');

ske_pred = [ske_pred_x1,ske_pred_y1];
                