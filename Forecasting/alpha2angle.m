function ang = alpha2angle(alpha, eigenWorms,eigen_num, vec_len)
% 
% 
% 
% 
% Kezhi Jan 2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

len_vec = length(vec_len);

ang = alpha * eigenWorms(1:eigen_num,:);

len_eigenWorms = size(eigenWorms,2);

    vec_vec_len = zeros(len_vec,1);
    vec_vec_len(1) = vec_len(1);
    for ii =2:length(vec_len);
        vec_vec_len(ii) = vec_vec_len(ii-1)+vec_len(ii);
    end
    
    [~, ind2] = histc(vec_vec_len/vec_vec_len(end), (1:len_eigenWorms)/len_eigenWorms);
    
    ind2(ind2<1)=1;
    centr_angle_aft = ang(ind2);
    if sum(abs(centr_angle_aft)>pi)>0
        larg = find(centr_angle_aft>pi);
        centr_angle_aft(larg) = -2*pi+centr_angle_aft(larg);
        smal = find(centr_angle_aft<-pi);
        centr_angle_aft(smal) = 2*pi+centr_angle_aft(smal);
    end
ang = centr_angle_aft;
