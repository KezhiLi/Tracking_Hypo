function X= resample_particles_Worm_body_hypo(XX, L_log,sub_num)

% Calculating Cumulative Distribution

L = exp(L_log - max(L_log));
Q = L / sum(L, 2);
R = cumsum(Q, 2);

% Generating Random Numbers

N = size(XX,1) * size(XX, 2);
T = rand(1, N);

% Resampling

[~, II] = histc(T, R);
I = II+1;
%tabulate(I)

% find the top N/sub_num entries with highest samples
I_unq = unique(I);
hist_I_unq = histc(I,I_unq);
% sort the histogram
[sorted_hist_I, index] = sort(hist_I_unq,'descend'); 

sorted_I = I_unq(index(1:round(N/sub_num)));
% only keep the entires with highest samples
mask_sort = zeros(1,length(Q));
mask_sort(sorted_I) = 1;

% %% this part can be neglected, and use sorted_I instead of I_new in the last sentence
% % Recalculate the probability and cum probability
% Q_new = Q.*mask_sort;
% Q_new = Q_new/sum(Q_new);
% R_new = cumsum(Q_new, 2);
% 
% % reduce the second sampling size to N/sub_num
% N_new = N/sub_num;
% T_new = rand(1, N_new);
% 
% % Resampling
% 
% [~, I_new] = histc(T_new, R_new);
%%
% X = XX(:, I_new+1);

for ii = 1:size(XX,1) ;
X{ii} = XX{floor((sorted_I(ii)-1)/size(XX,2))+1,mod(sorted_I(ii)-1,size(XX,2))+1};
end

