function [vec_len_pred_adjust, ave_vec_len] = bound_vec(vec_len_pred_now, len_min, len_max)
% function to make sure that the total length of vectors/segments are
% within the boundary of [minimal length, maximal length]. Otherwise,
% increasing/decreasing the length within the bound, and the variation is
% added to every vector/segment averagely. 
%
% Input: 
% vec_len_pred_now: an M * 1 vector, is the length of every
%                           vector/segments
% len_min: a scalar, the minimal bound of total length of the worm   
% len_max: a scalar, the maximum bound of total length of the worm
%
% Output:  
% vec_len_pred_adjust: an M * 1 vector, is the length of every
%                       vector/segment after adjustment
% sum_vec_len: the total length of the vectors/segments, which reprents the
%               length of the worm
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 31/03/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% sum lengths prediction
sum_len_pred = sum(vec_len_pred_now);

% The total length change
%         len_short = sum_len(ii) - sum_len_pred;
len_short = 0;

% set a threshold of the skeleton lengths
if sum_len_pred <len_min
    len_short = len_short + (len_min - sum_len_pred);
elseif sum_len_pred >len_max
    len_short = len_short + (len_max - sum_len_pred);
end

% Compensate the length change to the last 3 segments near tail
vec_len_pred_now = vec_len_pred_now + len_short/length(vec_len_pred_now);
%         vec_len_pred_adjust = 0.8*vec_len_pred +0.2*ave_vec_len_pred;
vec_len_pred_adjust = vec_len_pred_now;

% indication use, show the total length after adjustment
sum_vec_len = sum(vec_len_pred_adjust);

% average vector length
ave_vec_len = sum_vec_len/length(vec_len_pred_adjust);
