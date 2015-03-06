function flag  = sort_shape(flag)
% Because the frenet frame vector N will change directions frequently, in
% order to sort the points on the boundary, we need to keep vector N
% pointing towards one side of the worm permanatly. To achieve this, once
% we find two vector N in two succesive skeleton points with opposite directions, we
% multiply -1 to the second vector N to make it point to the right
% direction. In this function, we achieve this by count the number of -1 in
% the flag. -1 means that we need to change the +- of all following points
% on the skeleton. This operation is done in the loop.
%
% Copyright: Kezhi Li, 09/12/2014 CSC, MRC, London
% 11/12/2014
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

flag_1 = flag<-0.5;
flag_init = ones(length(flag),1);

summ = 0;
%
for ii = 1: length(flag);
    summ = summ + flag_1(ii);
    flag_init(ii) = (-1)^(mod(summ,2));
    %flag_init(ii) = (-1)^(mod(sum(flag_1(1:ii)),2));
end
flag = flag_init;