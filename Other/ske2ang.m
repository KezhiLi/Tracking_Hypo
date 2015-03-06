function [pt_ske, vec_len, angle_ske] = ske2ang(ske, seg_len)
% Calclate the change angle based on the skeleton
% Input: ske: an M * 2 matrix, the skeleton (usually obtained froma spline)
%        seg:len : a scalar, which determines how many points in each
%        segment
% Output: pt_ske: an m * 2 matrix, discreted points on the skeleton, start from head, with head segment length equal or smaller than the seg_len 
%         vec_len: a scalar, the length of segement, start from head
%         angle_ske: a m * 1 vector, the angles between two adjacent segment/vector, 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 11/12/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


% ske2point, start from the tail
if seg_len > 1
    if mod(length(ske)-2,seg_len)>0
        pt_ske = ske([end:-seg_len:2, 1]',:);
    else
        pt_ske = ske([end:-seg_len:3, 1]',:);
    end
else
    pt_ske = ske([end:-1: 1]',:);
end
% change direction to starting from head
pt_ske = pt_ske(end:-1:1,:);

% calculate the segment vector between each two adjacent points
vector = pt_ske(1:end-1,:)-pt_ske(2:end,:);
% convension from cart 2 polar coordinates
[angle_ske, vec_len] = cart2pol(vector(:,1),vector(:,2));

% debug use
if length(angle_ske)==0
    angle_ske
end



% % calculate the real length of each vector
% vec_len = sqrt(vector(:,1).^2+vector(:,2).^2);
% 
% % vec a * vec b
% mul_vec = diag(vector(1:end-1,:)*vector(2:end,:)');
% % |vec a| * | vec b|
% abs_mul_vec = vec_len(1:end-1,:).*vec_len(2:end,:);
% 
% % vec a * vec b / (|vec a| * | vec b|)
% cos_alpha = mul_vec./abs_mul_vec; 
% % calculate the angle
% angle_ske(:,1) = acos(cos_alpha);
% angle_ske(:,2) = sign(mul_vec);






