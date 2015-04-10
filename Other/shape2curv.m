function [worm_shape_x, worm_shape_y] = shape2curv(worm_shape1, para_spl, t,ts, siz)
% shape to curve
% Input: worm_shape1: a ?*2 matrix, points on the shape (? is the length of points on shape)
%        para_spl: a scalar between [0,1], the spline parameter, a tradeoff between spline and
%                  linear curves 
%        t: a vector,
%        ts: a longer vector than 't', with points located nearly
%        siz: size of the image
% Output: [worm_shape_x, worm_shape_y]: the x indexes and y indexes of
%                                       points on the curve
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 29/01/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


if para_spl == 1
    worm_shape_xy =  round(interp1(t,worm_shape1,ts,'pchip'));  % spline, pchip(68s)
else
para_lin = 1 - para_spl;
% Calculate the curve
worm_shape_xy =  round(para_spl*interp1(t,worm_shape1,ts,'spline')+para_lin*interp1(t,worm_shape1,ts,'linear') );
end

% Set threshold
worm_shape_xy(worm_shape_xy<1)=1;

% Set boundary inside the image
worm_shape_xy(worm_shape_xy(:,1)>siz(2),1)=siz(2);
worm_shape_xy(worm_shape_xy(:,2)>siz(1),2)=siz(1);

% Output values
worm_shape_xy = round(worm_shape_xy);

worm_shape_x = worm_shape_xy(:,1);
worm_shape_y = worm_shape_xy(:,2);

% % delete duplicated continues points
% ind_haf = round(size(worm_shape_xy,1)/2);
% worm_shape_xy_unq1 = unique(worm_shape_xy(1:ind_haf,:),'rows','stable');
% worm_shape_xy_unq2 = unique(worm_shape_xy(ind_haf+1:end,:),'rows','stable');
% worm_shape_xy_unq = [worm_shape_xy_unq1;worm_shape_xy_unq2];
% 
% worm_shape_x = worm_shape_xy_unq(:,1);
% worm_shape_y = worm_shape_xy_unq(:,2);




        
        