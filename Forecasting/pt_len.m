function [len, dis] = pt_len(Curv)
% 
% 
% 
% 
% 

pt_st = Curv(1:end-1,:);
pt_ed = Curv(2:end,:);

len = pt_ed - pt_st;

if size(len)==[0,0]
    Curv
end

dis = sqrt(len(:,1).^2+len(:,2).^2);
len = sum(dis);



