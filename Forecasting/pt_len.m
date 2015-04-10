function [len, dis] = pt_len(Curv)
% calculate the total length according the points on the curve
% input:  Curv: an M*2 vector, is the set of points on the Curv, from one
%               end to another end
% Output: len: the total length of the curv
%         dis: the vector, represents the length of segment between points
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 10/04/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


pt_st = Curv(1:end-1,:);
pt_ed = Curv(2:end,:);

len = pt_ed - pt_st;

if size(len)==[0,0]
    Curv
end

dis = sqrt(len(:,1).^2+len(:,2).^2);
len = sum(dis);



