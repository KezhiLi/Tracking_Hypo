function texture_output = texture_bef(k_5_fold, num1, mid_num, X_old_1_xy, width, Y_k_gray, texture)
% 
% 
% 
% 
% 
% 

shrink_para = 0.8;

    texture_newY{k_5_fold} = ske2tex(X_old_1_xy, width*shrink_para, Y_k_gray);
    num_pt_text = size(texture_newY{k_5_fold},1);
    texture_mtx(k_5_fold,1:2)=sum(abs(texture_newY{k_5_fold}-texture{k_5_fold-1}))/num_pt_text;
   % text_mtx_now = texture_mtx(k_5_fold,1:2)

    seq_cor1 = crossCheck(texture_newY{k_5_fold}(:,1), texture{k_5_fold-1}(:,1), num1);
    texture_mtx(k_5_fold,3) = find( min(seq_cor1)==seq_cor1)-mid_num;

texture_output{1} = texture_newY{k_5_fold};
texture_output{2} = texture_mtx(k_5_fold,1:2);