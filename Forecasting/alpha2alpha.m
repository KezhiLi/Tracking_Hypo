function [alpha_new] = alpha2alpha(alpha, jj)
% 
% 
% 
% 
% 
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

alpha_new = alpha - 0.2;
chg = str2num(dec2base(jj-1,3));

% if length(alpha)==4;
%     chg = chg - '1111';
% elseif length(alpha)==5;
%     chg = chg - '11111';
% elseif length(alpha)==6;
%     chg = chg - '111111';
% else
%     print('error alpha2alpha')
% end
if jj>1
    for ii = 1:length(chg);
        alpha_new(end+1-ii) = alpha_new(end+1-ii) + chg(end+1-ii)*0.2;
    end
end






