function [projectedAmps, fit_cuv_NonSaLeng] = eigenWormProject_NonSaLeng(eigenWorms, angleArray, numEigWorms,  varargin)
% eigenworm projection, for two vectors with different lengths
% 
% 
% 
% 
%     
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 15/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

if nargin >3
    vec_len = varargin{1}; 
    vec_vec_len = zeros(size(vec_len));
    vec_vec_len(1) = vec_len(1);
    for ii =2:length(vec_len);
        vec_vec_len(ii) = vec_vec_len(ii-1)+vec_len(ii);
    end
end

    centr_ang_chg = angleArray(1:end-1)-angleArray(2:end);
    % check if we have sudden change in angles
    coord_chg = sum(abs(centr_ang_chg)> 4);
    coord_chg_ind = find(abs(centr_ang_chg)> 4); 
    centr_angle_conti = angleArray;
    if coord_chg>0;
        for n_chg = 1:coord_chg; 
%             debug use only
%             if length(centr_ang_chg) == coord_chg_ind(n_chg)
%                 length(centr_ang_chg)
%                 coord_chg_ind(n_chg)
%             end
            centr_angle_conti(coord_chg_ind(n_chg)+1:end)=angleArray(coord_chg_ind(n_chg)+1:end)....
                +2*pi*(sign(angleArray(coord_chg_ind(n_chg))-angleArray(coord_chg_ind(n_chg)+1)));
       end
    end


    mean_centr_angle = mean(centr_angle_conti);
    centr_angle_m = centr_angle_conti - mean_centr_angle;
    
%project a given angle array onto predefined eigenWorms.  Only use the
%first numEigWorms eigenWorms.
len_eigenWorms = size(eigenWorms,2);
len_angleArray = size(centr_angle_m, 1);
num_angleArray = size(centr_angle_m, 2);
projectedAmps = NaN(num_angleArray, numEigWorms);


%Calculate time series of projections onto eigenworms
for i = 1:num_angleArray;
    if nargin ==3;
        ind = round((1:len_eigenWorms)/len_eigenWorms*len_angleArray);
    else
        [~, ind] = histc((1:len_eigenWorms)/len_eigenWorms, vec_vec_len/vec_vec_len(end));
    end
    ind(ind<1)=1;
    rawAngles = centr_angle_m(ind,i);
    for j = 1:numEigWorms
        projectedAmps(i,j) = eigenWorms(j,:)*rawAngles;
    end
end

    centr_angle_adj_mean = projectedAmps*eigenWorms(1:numEigWorms,:)+ mean_centr_angle;
    
    if nargin ==3;
        ind2 = round((1:len_angleArray)/len_angleArray*len_eigenWorms);
    else
        [~, ind2] = histc(vec_vec_len/vec_vec_len(end), (1:len_eigenWorms)/len_eigenWorms);
    end
    ind2(ind2<1)=1;
    centr_angle_aft = centr_angle_adj_mean(ind2);
    if sum(abs(centr_angle_aft)>pi)>0
        larg = find(centr_angle_aft>pi);
        centr_angle_aft(larg) = -2*pi+centr_angle_aft(larg);
        smal = find(centr_angle_aft<-pi);
        centr_angle_aft(smal) = 2*pi+centr_angle_aft(smal);
    end
    
fit_cuv_NonSaLeng = centr_angle_aft;