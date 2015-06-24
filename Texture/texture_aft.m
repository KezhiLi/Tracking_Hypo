function texture_output2 = texture_aft(k_5_fold, X_1_xy, width, Y_k_gray, texture, texture_newY, num1, mid_num  )    
% 
% 
% 
% 
% 
% 
% 

shrink_para = 0.8;

texture{k_5_fold} = ske2tex(X_1_xy, width*shrink_para, Y_k_gray);
num_pt_text = size(texture{k_5_fold},1);
texture_mtx(k_5_fold,4:5)=sum(abs(texture{k_5_fold}-texture_newY{k_5_fold}))/num_pt_text;
    

seq_cor2 = crossCheck(texture{k_5_fold}(:,1), texture_newY{k_5_fold}(:,1), num1);
texture_mtx(k_5_fold,6) = find( min(seq_cor2)==seq_cor2)-mid_num;
    
texture_mtx(k_5_fold,7:8)=sum(abs(texture{k_5_fold}-texture{k_5_fold-1}))/num_pt_text;
seq_cor3 = crossCheck(texture{k_5_fold}(:,1), texture{k_5_fold-1}(:,1), num1);
texture_mtx(k_5_fold,9) = find( min(seq_cor3)==seq_cor3)-mid_num; 

texture_output2{1} = texture;
texture_output2{2} = texture_mtx;

    
    

