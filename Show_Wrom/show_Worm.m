function show_Worm(worm_show, Y_k, width, seg_len, ind, diff, loc, inn_result)
% show the worm in image
% Input: worm_show: a ?*2 vector, the points on worm's skeleton
%        Y_k: a matrix, current real image
%        width: a scalar, the width(pixels) of the worm (radius)
%        seg_len: a scalar, the length of segment
% Output: an image(figure) shown in matlab, can be saved later
% 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 23/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%% Show current frame of video
figure(1)
if loc == 1;
image(Y_k)
%imshow(Y_k)
title('+++ Tracking Single Worm +++')
end

c=['diff ',num2str(loc), ' is ',num2str(round(diff))];
text(5,2+loc*10,c);

if ind == 1
    c_inn_result = ['average diff 1 is ',num2str(mean(inn_result(inn_result>0)))];
    text(5,2+(loc+1)*10,c_inn_result);
end

%% Calculate the skeleton

% the length of points on skeleton
m_fre_pt = size(worm_show,1);
% vector prepared for skeleton
t = 1:m_fre_pt ;
ts = 1:1/(seg_len*2):m_fre_pt;

worm_ske1 = worm_show;
% The skeleton curve, which is a combination of spline curve and linear curve
worm_ske_curv_x =  0.7*interp1(t,worm_ske1(:,1),ts,'spline')+0.3*interp1(t,worm_ske1(:,1),ts,'linear') ;
worm_ske_curv_y =  0.7*interp1(t,worm_ske1(:,2),ts,'spline')+0.3*interp1(t,worm_ske1(:,2),ts,'linear') ;

%% Based on the skeleton curve and width, calculate the contour on both sides
% of the worm. 

% % vector prepared for contour   
% tc = 1:2*m_fre_pt ;
% tsc = 1:1/(seg_len*2):2*m_fre_pt;
% 
% % Calculate the perpendicular direction of each points on skeleton
% [worm_ske1_T,worm_ske1_N] = frenet_TN(worm_ske1(:,1),worm_ske1(:,2));
% % Calculate points on the contour
% [worm_shape,worm_body] = ske2shape(worm_ske1, worm_ske1_N, width, -0.4);
% %worm_shape1 = round(worm_shape);
% worm_shape1 = worm_shape;
% % Calculate the curve links all the points on contour
% worm_shape_x =  0.7*interp1(tc,worm_shape1(:,1),tsc,'spline')+0.3*interp1(tc,worm_shape1(:,1),tsc,'linear') ;
% worm_shape_y =  0.7*interp1(tc,worm_shape1(:,2),tsc,'spline')+0.3*interp1(tc,worm_shape1(:,2),tsc,'linear') ;

%% Show the skeleton and contour on the current figure

hold on
%plot(worm_shape_x,worm_shape_y,'LineWidth',2,'Color',[1 0.1*ind 0.16*ind]);
plot(worm_ske_curv_x(5:end-4),worm_ske_curv_y(5:end-4),'LineWidth',2,'Color',[1 0.1*ind 0.16*ind]);

hold off 

