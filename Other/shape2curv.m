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

para_lin = 1 - para_spl;

% Calculate the curve
worm_shape_x =  round(para_spl*interp1(t,worm_shape1(:,1),ts,'spline')+para_lin*interp1(t,worm_shape1(:,1),ts,'linear') );
worm_shape_y =  round(para_spl*interp1(t,worm_shape1(:,2),ts,'spline')+para_lin*interp1(t,worm_shape1(:,2),ts,'linear') );

% Set threshold
worm_shape_x(worm_shape_x<1)=1;
worm_shape_y(worm_shape_y<1)=1;

% Set boundary inside the image
worm_shape_x(worm_shape_x>siz(2))=siz(2);
worm_shape_y(worm_shape_y>siz(1))=siz(1);

% Output integer values
worm_shape_x = round(worm_shape_x);
worm_shape_y = round(worm_shape_y);
        
        