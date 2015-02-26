clc;
clear;

rgb = imread('P:\Kezhi\FromAndre\images\DCR_F_618_13065_3.bmp');
%rgb = imread('P:\Kezhi\MyCode\Sample\video\sample1\o_1235793177532ed3-050.jpg');
%rgb = imread('C_elegan2.jpg');
I = rgb2gray(rgb);
%I = rgb;
[mm,nn]= size(I);
%figure, imshow(I), title('Original')
II = 255 - I;

II(II<10)=(median(im2double(reshape(II,size(II,1)*size(II,2),1)))*255);

se = strel('disk', 1);
II = imclose(II, se);
%figure
%imshow(II), title('Opening-closing (Ioc II)')

 III = imclose(II, se);
% figure
% imshow(III), title('Opening-closing (Ioc III)')

%figure, hist(double(I),256);
level = graythresh(II);
BW = im2bw(II,level);
% figure, imshow(BW)

reBW = BW;
% reBW = -im2bw(I,level)+1;
%  figure, imshow(reBW)

log_reBW = logical(reBW);

BWdfill = imfill(log_reBW, 'holes');
% figure, imshow(BWdfill);
% title('binary image with filled holes');

BW2 = bwareaopen(BWdfill, 50); 
% figure, imshow(BW2);
% title('remove small areas');


BWoutline = bwperim(BW2);
figure, imshow(BWoutline);

%[NN, curv, Verticles, Lines] = curvature(BWoutline);
[NN, curv, Vertices, Lines, Vertices_old] = curvature_areainput(BW2);

[mm_Ver,nn_Ver] = size(Vertices); 
curv_norm=curv(:,1).^2+curv(:,2).^2;
curv_point1 =  find(curv_norm(:)==max(curv_norm));
curv_point1 = curv_point1(1);
curv_norm_double=[curv_norm;curv_norm];
curv_norm_opp = curv_norm_double(curv_point1+round(mm_Ver*0.2):curv_point1+round(mm_Ver*0.8));   % search the opposite side of the contour
curv_point2 =  find(curv_norm(:)==max(curv_norm_opp));
curv_point2 = curv_point2(1);
figure, plot(curv_norm)

Segout = gray2rgb(I);
Segout(BWoutline) = 1;

V_x=round([Vertices(curv_point1,1),Vertices(curv_point2,1)]);
V_y=mm-round([Vertices(curv_point1,2),Vertices(curv_point2,2)]);

for ii = V_y(1)-1-2:V_y(1)-1+2;
    for jj = V_x(1)-1-2:V_x(1)-1+2;
        Segout(ii,jj,:) = [0.2,0.2,1];
    end
end
for ii = V_y(2)-1-2:V_y(2)-1+2;
    for jj = V_x(2)-1-2:V_x(2)-1+2;
        Segout(ii,jj,:) = [0.2,0.2,1];
    end
end

% debug use
%hold on, plot([Verticles(curv_point1,1),Verticles(curv_point2,1)],[mm-Verticles(curv_point1,2),mm-Verticles(curv_point2,2)],'bo','MarkerSize', 10);

% The worm is roughly divided into 24 segments of musculature (i.e., hinges
% that represent degrees of freedom) on each side. Therefore, 48 segments
% around a 2-D contour.
sWormSegs = 24;
cWormSegs = 2 * sWormSegs;

% Clean up the worm's contour.
contour = cleanWorm(Vertices, size(Vertices, 1) / cWormSegs);
% contour = Vertices;

% Compute the contour's local high/low-frequency curvature.
% Note: worm body muscles are arranged and innervated as staggered pairs.
% Therefore, 2 segments have one theoretical degree of freedom (i.e. one
% approximation of a hinge). In the head, muscles are innervated
% individually. Therefore, we sample the worm head's curvature at twice the
% frequency of its body.
% Note 2: we ignore Nyquist sampling theorem (sampling at twice the
% frequency) since the worm's cuticle constrains its mobility and practical
% degrees of freedom.
cCCLengths = circComputeChainCodeLengths(contour);
wormSegLength = (cCCLengths(1) + cCCLengths(end)) / cWormSegs;
hfAngleEdgeLength = wormSegLength;
hfCAngles = circCurvature(contour, hfAngleEdgeLength, cCCLengths);
lfAngleEdgeLength = 2 * hfAngleEdgeLength;
lfCAngles = circCurvature(contour, lfAngleEdgeLength, cCCLengths);

% Compute the contour's local high/low-frequency curvature maxima.
%[mhfCMaxP mhfCMaxI] = maxPeaksCircDist(mhfCAngles, hfAngleEdgeLength, ...
%    cCCLengths);
[lfCMaxP lfCMaxI] = maxPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
    cCCLengths);


 % Compute the contour's local low-frequency curvature minima.
 [lfCMinP lfCMinI] = minPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
        cCCLengths);
 
 headI = curv_point1;
 tailI = curv_point2;
 % Compute the worm's skeleton.
 [skeleton cWidths] = linearSkeleton(headI, tailI, lfCMinP, lfCMinI, ...
        lfCMaxP, lfCMaxI, contour, wormSegLength, cCCLengths);


    % Measure the skeleton's chain code length.
    sCCLengths = computeChainCodeLengths(skeleton);
    sLength = sCCLengths(end);

    % Measure the skeleton angles (curvature).
    lfAngleEdgeLength = sCCLengths(end) * 2 / sWormSegs;
    sAngles = curvature(skeleton, lfAngleEdgeLength, sCCLengths);

 skeleton_sw = [mm-skeleton(:,2),skeleton(:,1)];
 Img_ske = zeros(mm,nn);
 skeleton_sw_1d = sub2ind(size(Img_ske),  skeleton_sw(:,1),  skeleton_sw(:,2));
 
  % reverse skeleton to skeleton_sw to make it good for plot
   skeleton_sw = [mm-skeleton(:,2),skeleton(:,1)];
   % find the central point of the skeleton
   hf_ske_sw = skeleton_sw(round(length(skeleton)/2),:);
   
   % To show the skeleton
   Img_ske = zeros(mm,nn);
   skeleton_sw_1d = sub2ind(size(Img_ske),  skeleton_sw(:,1),  skeleton_sw(:,2));
    % as the red colour
   Img_ske(skeleton_sw_1d) = 1;
   Segout(logical(Img_ske)) = 0.8;
 
   % draw the central point of skeleton as a green point
   Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,2) = 1;
   Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,1) = 0.2;
   Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,3) = 0.2;

 
 
figure, imshow(Segout), title('outlined original image');











