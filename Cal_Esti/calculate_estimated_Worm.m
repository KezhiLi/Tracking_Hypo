function [worm_show,X,inn_result_k] = calculate_estimated_Worm(X, worm_sav,  C_k, width, tt, seg_len, ind, size_blk)
% function to calculate the estimated worm
% Input:    X: a cell matrix, all good hypotheses
%           worm_sav: a ?*2 vector, the skeleton saved before
%           C_k: a matrix, current frame
%           width: a scalar, the worm's width
%           tt: a scalar, the averaging number
%           seg_len: a scalar, the length(points) of segment
% Output:   worm_show: a ?*2 vector, the skeleton will be shown
%           X: a cell matrix, new generated with worm_shown in it
%           inn_result_k: the coherence of the worm_show and C_k(after thresholding)
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 24/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% height and width
Npix_h = size(C_k, 1);
Npix_w = size(C_k, 2);

mask_ori = ones(Npix_h,1)*ones(1,Npix_w);

% empty image, faster if they are generated in these ways
area_hypo = zeros(Npix_h,0)*zeros(0,Npix_w);
sav_area_hypo = zeros(Npix_h,0)*zeros(0,Npix_w);

% length of skeleton points
m_fre_pt = size(X{1,1}.xy,1);
% points for extension in contour
t = 1:2*m_fre_pt ;
ts = 1:1/(seg_len*2):2*m_fre_pt;
% points for extension in body, which consist of total_seg (eg. 9) curves on and inside
% the body contour
div = floor(width + 1);
total_seg = (div+1)*2-1;
t2 = 1:total_seg*m_fre_pt ;
ts2 = 1:1/(seg_len*2):total_seg*m_fre_pt;

worm_ske = X{ind}.xy;
% result is calculated from weighted averaging top tt hypotheses 
para = exp(-(1:tt)*0.2+0.2);   % weights are exponentially decaying
sum_para = sum(para);
for ii=2:tt;
    worm_ske = worm_ske + X{ii}.xy*para(ii);
end
%worm_ske1 = round(worm_ske/sum_para);
worm_ske1 = worm_ske/sum_para;

%worm_ske1 = worm_ske;
sum_C_k = sum(sum(C_k))/700;

%% Calculate the new worm

% calculate frenet vectors
[worm_ske1_T,worm_ske1_N] = frenet_TN(worm_ske1(:,1),worm_ske1(:,2));
% calculate points on contour
[worm_shape1,worm_body] = ske2shape(worm_ske1, worm_ske1_N, width, -0.40);
%[worm_shape_x, worm_shape_y] = shape2curv(worm_shape1, 0.5, t, ts, size(area_hypo));

% calculate points inside body
[worm_body_x, worm_body_y] = shape2curv(worm_body, 1, t2, ts2, size(area_hypo));   % 0.7

% debug use
if isnan(sum(sum(worm_body_x)))
    worm_body_x
end
       
%area_hypo_1d = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );

mask_head = square_mtx_fast(mask_ori, X{ii}.xy(end-1,:), size_blk);
mask = square_mtx_fast(mask_head, X{ii}.xy(2,:), size_blk);    
                
area_hypo_1d_inside = sub2ind(size(area_hypo),   worm_body_y , worm_body_x );
        
I = (min(worm_body_y) >= 1 & max(worm_body_y) <= Npix_h);
J = (min(worm_body_x) >= 1 & max(worm_body_x) <= Npix_w);
    
if I && J
        
    %debug use
    if (sum(abs(area_hypo_1d_inside-round(area_hypo_1d_inside)))>0)||(min(min(area_hypo_1d_inside))<1)
            area_hypo_1d_inside
    end
        
    area_hypo(area_hypo_1d_inside) = 1;

    area_hypo1 = area_hypo;
    log_reBW_hypo = logical(area_hypo1);
        
    BWdfill_hypo = log_reBW_hypo;

    % open operator
    BW2_hypo = bwareaopen(BWdfill_hypo, 50); 
        
        % calculate the coherence for current hypothese
    D_new = sum(sum(imabsdiff(C_k,BW2_hypo).*mask))/sum_C_k;           

end
%% Calculate the saved worm

[worm_sav_ske1_T,worm_sav_ske1_N] = frenet_TN(worm_sav(:,1),worm_sav(:,2));
[worm_sav_shape1,worm_sav_body] = ske2shape(worm_sav, worm_sav_ske1_N, width, -0.40);

%[worm_sav_shape_x, worm_sav_shape_y] = shape2curv(worm_sav_shape1, 0.5, t, ts, size(sav_area_hypo));

[worm_sav_body_x, worm_sav_body_y] = shape2curv(worm_sav_body, 1, t2, ts2, size(sav_area_hypo)); % 0.7
         
%sav_area_hypo_1d = sub2ind(size(sav_area_hypo),   worm_sav_shape_y , worm_sav_shape_x );
                
sav_area_hypo_1d_inside = sub2ind(size(sav_area_hypo),   worm_sav_body_y , worm_sav_body_x );
        
II = (min(worm_sav_body_y) >= 1 & max(worm_sav_body_y) <= Npix_h);
JJ = (min(worm_sav_body_x) >= 1 & max(worm_sav_body_x) <= Npix_w);
    
if II && JJ
        
    %sav_area_hypo(sav_area_hypo_1d) = 1;
        
    sav_area_hypo(sav_area_hypo_1d_inside) = 1;

    sav_area_hypo1 = sav_area_hypo;
    sav_log_reBW_hypo = logical(sav_area_hypo1);
        
    sav_BWdfill_hypo = sav_log_reBW_hypo;

    sav_BW2_hypo = bwareaopen(sav_BWdfill_hypo, 50); 
    
    % Calculate the coherence of the old worm
    D_old = sum(sum(imabsdiff(C_k,sav_BW2_hypo).*mask))/sum_C_k;           
end
    
    % if the new worm is better than the old worm
    if D_new<D_old
        % show a combination of the new worm and old worm
        worm_show = round(0.8*worm_ske1+0.2*worm_sav);
        % save the coherence of the new worm
        inn_result_k = D_new;
        % save new worm to X
        X{end}.xy = worm_ske1;
        % save the shown worm to X
        X{end-1}.xy = worm_show;
    else % if the new worm is worse than the old worm
        % show the old worm
        worm_show = worm_sav;
        % save the coherence of the old worm
        inn_result_k = D_old;
        % save the old worm to X
        X{end}.xy = worm_show;
        % save the new worm to X
        X{end-1}.xy = worm_ske1;
    end

