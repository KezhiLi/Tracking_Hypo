function [worm_show,X,inn_result_k] = calculate_estimated_Worm_body_hypo(X, worm_sav,  C_k, width, tt)

area_hypo = zeros(size(C_k,1),0)*zeros(0,size(C_k,2));
sav_area_hypo = zeros(size(C_k,1),0)*zeros(0,size(C_k,2));

Npix_h = size(C_k, 1);
Npix_w = size(C_k, 2);

%X_mean = mean(X, 2);

    seg_len = 8;
    m_fre_pt = size(X{1,1}.xy,1);
    t = 1:2*m_fre_pt ;
    ts = 1:1/(seg_len*2):2*m_fre_pt;
    t2 = 1:9*m_fre_pt ;
    ts2 = 1:1/(seg_len*2):9*m_fre_pt;

worm_ske = X{1}.xy;
% result is calculated from weighted averaging top tt hypotheses 
%tt=5;
para = exp(-(1:tt)*0.2+0.2);
sum_para = sum(para);
for ii=2:tt;
    worm_ske = worm_ske + X{ii}.xy*para(ii);
end
worm_ske1 = round(worm_ske/sum_para);

sum_C_k = sum(sum(C_k));
%%

[worm_ske1_T,worm_ske1_N] = frenet_TN(worm_ske1(:,1),worm_ske1(:,2));
[worm_shape1,worm_body] = ske2shape(worm_ske1, worm_ske1_N, width, -0.40);

%[worm_shape_x, worm_shape_y] = shape2curv(worm_shape1, 0.5, t, ts, size(area_hypo));

[worm_body_x, worm_body_y] = shape2curv(worm_body, 0.7, t2, ts2, size(area_hypo));
if isnan(sum(sum(worm_body_x)))
    worm_body_x
end
       
       
%area_hypo_1d = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );
                
area_hypo_1d_inside = sub2ind(size(area_hypo),   worm_body_y , worm_body_x );
        
    I = (min(worm_body_y) >= 1 & max(worm_body_y) <= Npix_h);
    J = (min(worm_body_x) >= 1 & max(worm_body_x) <= Npix_w);
    
    if I && J
        
        %area_hypo(area_hypo_1d) = 1;
        if (sum(abs(area_hypo_1d_inside-round(area_hypo_1d_inside)))>0)||(min(min(area_hypo_1d_inside))<1)
            area_hypo_1d_inside
        end
        area_hypo(area_hypo_1d_inside) = 1;

        area_hypo1 = area_hypo;
        log_reBW_hypo = logical(area_hypo1);
        
        BWdfill_hypo = log_reBW_hypo;

        BW2_hypo = bwareaopen(BWdfill_hypo, 50); 

        D_new = sum(sum(imabsdiff(C_k,BW2_hypo)))/sum_C_k;           
        %D = sum(sum(abs(C-BW2_hypo)));

    end
%%

[worm_sav_ske1_T,worm_sav_ske1_N] = frenet_TN(worm_sav(:,1),worm_sav(:,2));
[worm_sav_shape1,worm_sav_body] = ske2shape(worm_sav, worm_sav_ske1_N, width, -0.40);

%[worm_sav_shape_x, worm_sav_shape_y] = shape2curv(worm_sav_shape1, 0.5, t, ts, size(sav_area_hypo));

[worm_sav_body_x, worm_sav_body_y] = shape2curv(worm_sav_body, 0.7, t2, ts2, size(sav_area_hypo));
         
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

        D_old = sum(sum(imabsdiff(C_k,sav_BW2_hypo)))/sum_C_k;           
        %D = sum(sum(abs(C-BW2_hypo)));

    end
    
    
    if D_new<D_old
        worm_show = round(0.7*worm_ske1+0.3*worm_sav);
        inn_result_k = D_new;
        X{end}.xy = worm_ske1;
        X{end-1}.xy = worm_show;
    else 
        worm_show = worm_sav;
        inn_result_k = D_old;
        X{end}.xy = worm_show;
        X{end-1}.xy = worm_ske1;
    end
            









%drawnow
