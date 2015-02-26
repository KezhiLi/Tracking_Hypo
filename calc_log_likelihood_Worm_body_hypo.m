function L = calc_log_likelihood_Worm_body_hypo(Xstd_rgb, Xrgb_trgt, XX, Y, width) %#codegen

Npix_h = size(Y, 1);
Npix_w = size(Y, 2);

Npop_particles = size(XX,1);
sub_num = size(XX,2);

N = Npop_particles*sub_num;

L = zeros(1,N);
%Y = permute(Y, [3 1 2]);
Y = rgb2gray(Y);

[mm,nn]= size(Y);
%figure, imshow(I), title('Original')
II = 255 - Y;

 se = strel('disk', 1);
% II = imclose(II, se);

%figure, hist(double(I),256); 

% we set an artificial parameter 0.92 here
level = graythresh(II)*0.90;
BW = im2bw(II,level);
%figure, imshow(BW)

reBW = BW;
% reBW = -im2bw(I,level)+1;
% figure, imshow(reBW)

log_reBW = logical(reBW);

BWdfill = log_reBW;
%BWdfill = imfill(log_reBW, 'holes');
%figure, imshow(BWdfill);
%title('binary image with filled holes');

C = bwareaopen(BWdfill, 50); 
%figure, imshow(BW2);
%title('remove small areas');



A = -log(sqrt(2 * pi) * Xstd_rgb);
B = - 0.5 / (Xstd_rgb.^2);

%XX = round(XX);
    
    seg_len = 8;
    m_fre_pt = size(XX{1,1}.xy,1);
    t = 1:2*m_fre_pt ;
    ts = 1:1/(seg_len*2):2*m_fre_pt;
    
    limit = 3;

for ii = 1:Npop_particles;
    for jj = 1:sub_num;
       
        area_hypo = zeros(size(Y,1),0)*zeros(0,size(Y,2));
        
        worm_contour1 = round(ske2shape(XX{ii,jj}.xy, XX{ii,jj}.N, width, -0.20));
        
        worm_shape_x =  round(0.5*interp1(t,worm_contour1(:,1),ts,'spline')+0.5*interp1(t,worm_contour1(:,1),ts,'linear') );
        worm_shape_y =  round(0.5*interp1(t,worm_contour1(:,2),ts,'spline')+0.5*interp1(t,worm_contour1(:,2),ts,'linear') );
        
        worm_shape_x(worm_shape_x<1)=1;
        worm_shape_y(worm_shape_y<1)=1;
         
        area_hypo_1d = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );
        
        for tt= 1:width-limit;
            worm_contour{tt}=round(ske2shape(XX{ii,jj}.xy, XX{ii,jj}.N, tt, -0.20));
            
            worm_shape_x =  round(0.5*interp1(t,worm_contour{tt}(:,1),ts,'spline')+0.5*interp1(t,worm_contour{tt}(:,1),ts,'linear') );
            worm_shape_y =  round(0.5*interp1(t,worm_contour{tt}(:,2),ts,'spline')+0.5*interp1(t,worm_contour{tt}(:,2),ts,'linear') );
        
            worm_shape_x(worm_shape_x<1)=1;
            worm_shape_y(worm_shape_y<1)=1;
            
            area_hypo_1d_inside{tt} = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );
        end
        
    
    I = (min(worm_shape_y) >= 1 & max(worm_shape_y) <= Npix_h);
    J = (min(worm_shape_x) >= 1 & max(worm_shape_x) <= Npix_w);
    
    if I && J
        
        area_hypo(area_hypo_1d) = 1;
        
        for tt=1:width-limit;
            area_hypo(area_hypo_1d_inside{tt}) = 1;
        end

         %area_hypo1 = imclose(area_hypo, se);
        area_hypo1 = area_hypo;
        log_reBW_hypo = logical(area_hypo1);
        
        % cancel it to save time
        % BWdfill_hypo = imfill(log_reBW_hypo,  fliplr(worm_contour1(round(size(worm_contour1,1)/2),:))  );
        BWdfill_hypo = log_reBW_hypo;
        
        %figure, imshow(BWdfill);
        %title('binary image with filled holes');

        BW2_hypo = bwareaopen(BWdfill_hypo, 50); 
        %figure, imshow(BW2);
        %title('remove small areas');

          D = sum(sum(imabsdiff(C,BW2_hypo)));           
         %D = sum(sum(abs(C-BW2_hypo)));

        D2 = D' * D;
        
        L((ii-1)*sub_num+jj) =  A + B * D2;

         %L(k) =  1/sqrt((good_pt(1)-m)^2+(good_pt(2)-n)^2); 
    else
        
        %L(k) = 0;
        L(k) = -Inf;
    end
    end
end
