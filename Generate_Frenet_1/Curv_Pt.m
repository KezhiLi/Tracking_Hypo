function [curv_norm, curv_point1, curv_point2] = Curv_Pt(Vertices, curv)
% CURV_PT computes the norm of curvature, and find the points of the
% maximum curvatures, which are deemed as head and tail
% 
% INPUT: Vertices: Vertices vector
%        curv:  ?*2 curv vector
% OUTPUT: curv_norm:  magnitudes(norm) of curvature
%        curv_point1, curv_point2: 1*2 vector, two points with maximum curvatures
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 27/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

[mm_Ver,nn_Ver] = size(Vertices); 
curv_norm=curv(:,1).^2+curv(:,2).^2;  % compute the magnitudes of curvature
curv_point1 =  find(curv_norm(:)==max(curv_norm)); % find the maximum curvature
curv_point1 = curv_point1(1);   % in case there are more than 1 maximum curvature magnitudes
curv_norm_double=[curv_norm;curv_norm]; % stack two curvature vectors
curv_norm_opp = curv_norm_double(curv_point1+round(mm_Ver*0.2):curv_point1+round(mm_Ver*0.8));   % search the opposite side of the contour
curv_point2 =  find(curv_norm(:)==max(curv_norm_opp)); % find another local maximum in the other side 
curv_point2 = curv_point2(1);   % save the other local maximum curvature