function [worm_shape, worm_body] = ske2shape(xy, NN, W_bar, para )
%%
% Calculate the worm shape given skeleton from the head point
% Input:  xy: M * 2 matrix
%         NN: M * 2(or 3) matrix
%         W_bar: the estimated widest width of the worm
%         para: parameter that in the estimation equation, normally 0.1~0.2
% Output: worm_shape: L * 2 matrix, describle the points on the boundary of the worm
%
% Test data example:
% load Frenet_Pt
% xy = Frenet_Pt{2}.xy;
% NN = Frenet_Pt{2}.N(:,1:2);
% para = -0.40; 
% W_bar = 10;
%
% Copyright: Kezhi Li, 09/12/2014, CSC, MRC, London
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%%
% NN only has two dimention in our case
NN = NN(:,1:2);

% m+1 points in all
m = size(NN, 1) - 1;

% m+1 points, from 0 to m, apply the estimate equation
for ii = 0:m
    W(ii+1) = W_bar*(1 - exp(para*(m-abs(2*ii-m))));
end

% debug use: 
%figure, plot(0:m, W,'r')

% debug use:
% if min(abs(NN(:,1:2)))==[0 0]
%     NN
% end

% norm of NN, add a small value to avoid inf
%norm_NN = (ones(size(NN, 1),1)./(1e-10+sqrt(NN(:,1).^2+NN(:,2).^2) ))' ;
norm_NN = (ones(size(NN, 1),1)./(sqrt(NN(:,1).^2+NN(:,2).^2) ))' ;
% length of NN 
length_NN = length(norm_NN);

% for debug 
NN_old = NN;
norm_NN_old = norm_NN;

% in case NN = 0
% run this loop until all inf disappear
for kk=1:length_NN-1;
    % if there still inf, set it as a value as the next point on the
    % skeleton. If this point is the last point on the skeleton, set it as
    % a value as the previous point on the skeleton. 
%if sum(isinf(norm_NN))>0
if max(norm_NN)>1e9
   %ind_inf = find(isinf(norm_NN) == 1);
   ind_inf = find(norm_NN > 1e9);
   for ii = 1: length(ind_inf);
       if ind_inf(ii) == length_NN
          NN(ind_inf(ii),:) = NN(ind_inf(ii)-1,:);
          norm_NN(ind_inf(ii)) = norm_NN(ind_inf(ii)-1); 
       elseif ind_inf(ii) == 1;
          NN(ind_inf(ii),:) = NN(ind_inf(ii)+1,:);
          norm_NN(ind_inf(ii)) = norm_NN(ind_inf(ii)+1); 
       else
          NN(ind_inf(ii),:) = (NN(ind_inf(ii)+1,:)+NN(ind_inf(ii)-1,:)+1e-9)/2;
          if sum(abs(NN(ind_inf(ii),:)))< 1e-9;
              NN(ind_inf(ii),:) = NN(ind_inf(ii)+1,:);
          end
          norm_NN(ind_inf(ii)) = 1/(sqrt(NN(ind_inf(ii),1).^2+NN(ind_inf(ii),2).^2) );
       end
   end
else
    break;
end
end

% Normalize NN   
NN_normalized = [NN(:,1).*norm_NN' ,NN(:,2).*norm_NN'];

% sort the points on the boundary
width = [NN_normalized(:,1).*W',NN_normalized(:,2).*W'];
width_2 = [width(end,:);width(1:end-1,:)];

% to see if vector N change the direction dramatically
flag = sign(sum(width.*width_2,2));
% function created to sort the points on the shape boundary
flag = sort_shape(flag);

% width = mag * angle
% width = width.*kron(flag, ones(1,2));
width = width.*[flag,flag];

% calculate the points on the boundary on both sides
xy_N = xy + width;
xy_flipN = xy - width;

% % debug use
% if sum(sum(isnan(xy_N)))>0.5
%     xy_N
% end

% the output is a circle 
worm_shape = [xy_N; flipud(xy_flipN)];

%figure, plot(worm_shape(:,1),worm_shape(:,2),'r');

% fill the body of worm with points. 
% 0,0.25,0.5,0.75,1, five curves to fill the body

worm_body =worm_shape;

div = floor(W_bar + 1);
div_rev = 1/div;

for jj = (1-div_rev):-div_rev:div_rev;
    worm_body = [worm_body; xy + jj*width];
    worm_body = [worm_body; flipud(xy - jj*width)];
end
worm_body = [worm_body; xy];





        
    
