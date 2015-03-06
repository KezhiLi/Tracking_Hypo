function ske_hypo = ang2ske(ang_hypo, vec_len_pred, mid_hypo, mid_start)
% calculate the hypothesis skeleton given hypothesis angle, predicted
% vectors, hypothesis middle point, and the middle point index. The angle
% skeleton and angle vector are from head to tail. eg. the first value in
% ang_hypo denotes the angle from point2 pointing to point1, and the
% vec_len_pred is the vector from point2 pointing to point1. 
%
% Input: angle_hypo: an m * 1 vector, representing the hypothesis angle
%        vec_len_pred: an m * 2 matrix, representing the vectors link the sample points on
%        the skeleton. 
%        mid_hypo: hypothesis middle point
%        mid_start: the index of the middle point in ske_pt, which is 1
%        more element than ang_hypo
%Output: ske_hypo: an M * 2 matrix, which is the hypothesis skeleton
%         
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 15/12/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% change coordinates, from polar to cart
[XX, YY] =  pol2cart(ang_hypo, vec_len_pred);

% find the middle point, all other points are calculated based on the
% middle point as the starting point
ske_pt_hypo = zeros(size(ang_hypo,1)+1, 2);
ske_pt_hypo(mid_start,:) = mid_hypo;

% ske_pt_hypo: an (m+1) * 2 matrix, , which is the hypothesis skeleton
%         sample points

% calculate locations of the points from middle to the head
for ii = mid_start-1:-1:1;
    ske_pt_hypo(ii,:) = ske_pt_hypo(ii+1,:)+[XX(ii),YY(ii)];
end

% calculate locations of the points from middle to the tail
for ii = mid_start+1:size(ske_pt_hypo,1);
    ske_pt_hypo(ii,:) = ske_pt_hypo(ii-1,:)-[XX(ii-1),YY(ii-1)];
end

ske_hypo = ske_pt_hypo;




