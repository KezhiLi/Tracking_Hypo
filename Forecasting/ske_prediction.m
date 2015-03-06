function ske_pred = ske_prediction(speed, for_curv, back_curv, head, tail)
% Input: speed: a scalar \in (-infinity, infinity)
%        
%        
%        
% Output: ske_prediction: an L * 2 matrix, describle the points on the
%        predicted skeleton. 
% 
% Copyright: author: Kezhi Li, CSC, MRC, Imperial College, London
% 11/12/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% head = back_curv(end,:);
% tail = for_curv(end,:);


if speed>0;
    % find the index of head point in forward curve
    [neighborIds_for neighborDistances] = kNearestNeighbors(for_curv, head, 1);
    head_ind = neighborIds_for(1);
    
    % find the index of tail point in forward curve
    [neighborIds_for neighborDistances] = kNearestNeighbors(for_curv, tail, 2);
    % incase there are two points with 0 distance at the tailpoint
    tail_ind = neighborIds_for(1+(neighborDistances(2)==0));
    
    % genrate the predicted forward path with reverse points 
    path_for = for_curv(head_ind:-1:1,:);
    % calculate the moving distance for each predicted point
    [arclen,seglen] = arclength(path_for(:,1),path_for(:,2));
    cum_len = cumsum(seglen);
    % find the point on the forward path that is coincident to the forward moving speed  
    [pred_err, pred_ind] = min(abs(cum_len - speed));
    % locate the predicted head index in forward curve
    pred_head_ind = head_ind-pred_ind+1;
    % locate the predicted tail index in the forward curve
    pred_tail_ind = pred_head_ind + size(for_curv(head_ind:tail_ind,:),1)-1;
    % calculate the predicted skeleton based on the spline
    ske_pred = for_curv(pred_head_ind:pred_tail_ind,:);
    
elseif speed<0;
        % find the index of head point in forward curve
    [neighborIds_for neighborDistances] = kNearestNeighbors(back_curv, head, 2);
    % incase there are two points with 0 distance at the headpoint
    head_ind = neighborIds_for(1+(neighborDistances(2)==0));
    
    % find the index of tail point in backward curve
    [neighborIds_back neighborDistances] = kNearestNeighbors(back_curv, tail, 1);
    tail_ind = neighborIds_back(1);
    % genrate the predicted backward path 
    path_back = back_curv(tail_ind:-1:1,:);
    % calculate the moving distance for each predicted point
    [arclen,seglen] = arclength(path_back(:,1),path_back(:,2));
    cum_len = cumsum(seglen);
    % find the point on the backward path that is coincident to the backward moving speed  
    [pred_err, pred_ind] = min(abs(cum_len + speed));
    % locate the predicted tail index in backward curve
    pred_tail_ind = tail_ind-pred_ind+1;
    % locate the predicted head index in the backward curve
    pred_head_ind = pred_tail_ind + size(back_curv(tail_ind:head_ind,:),1)-1;
    % calculate the predicted skeleton based on the spline
    ske_pred = back_curv(pred_head_ind:-1:pred_tail_ind,:);
else
    % in the case that speed = 0, skeleton stays the same
    
    % find the index of head point in forward curve
    [neighborIds_for neighborDistances] = kNearestNeighbors(for_curv, head, 1);
    head_ind = neighborIds_for(1);
    
        % find the index of tail point in forward curve
    [neighborIds_for neighborDistances] = kNearestNeighbors(for_curv, tail, 1);
    tail_ind = neighborIds_for(1);
    
    ske_pred = for_curv(head_ind:tail_ind,:);
    
end