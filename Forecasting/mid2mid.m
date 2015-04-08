function [mid_hypo] = mid2mid(ske, pt_ske, mid_start, var_ske_direc, var_direc)
%  function used to calculate hypothesis middle point location given the
%  predicted middle point. It is obtained by adding a random index along the spline
%  direction and a random noise 
%
% Input: ske: 
%        an M * 2 matrix, representing the skeleton points location
%
%        pt_ske: 
%        an m * 2 matrix, representing the sample points on the ske
%        with segment length seg_len, the head of ske usually has a small
%        segement
%
%        mid_start: 
%        a scalar, the index of middle point in pt_ske
%
%        var_ske_direc, var_direc: 
%        the hypothesis is obtained from pridicted middle point plus a random index calculated from
%        var_ske_direc, and a random noise calculated from var_direc. 
%
% Output: mid_hypo:
%         a 1 * 2 vector, denoting the location of middle point 
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 15/12/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


pt_mid = pt_ske(mid_start,:);
pt_mid_mtx = [ones(size(ske,1),1)*pt_mid(1), ones(size(ske,1),1)*pt_mid(2)];
diff = sum(abs(ske-pt_mid_mtx),2);
[diff_seq diff_idx] = sort(diff,'ascend');

mid_ind = diff_idx(1);

% debug use
if length(mid_ind)>1
   mid_ind= mid_ind(1+(abs(mid_ind(1)-length(ske)/2)>abs(mid_ind(2)-length(ske)/2)));
end

mid_start_pred = ske(mid_ind+round(var_ske_direc*randn(1)),:);
mid_hypo = mid_start_pred + (var_direc * randn(1,2));


