function X = resample_particles_Worm(XX, L_log, fold)
% Function to resample the particles based on likelihood, and only save the
% particles with highest samples as X for the next iteration.
% 
% Input: XX: a cell matrix, all hypotheses 
%        L_log: the likelihood
% Output:X: a cell matrix, selected from XX that received the highest numbers
%           of samples 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 24/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

[N_particles, sub_num] = size(XX);

N_particles = N_particles * fold;

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

if length(index)< N_particles
    index = [index, ones(1,N_particles-length(index))*index(1)];
end
sorted_I = I_unq(index(1:N_particles));
% only keep the entires with highest samples
mask_sort = zeros(1,length(Q));
mask_sort(sorted_I) = 1;

for ii = 1:size(XX,1)*fold ;
X{ii} = XX{floor((sorted_I(ii)-1)/size(XX,2))+1,mod(sorted_I(ii)-1,size(XX,2))+1};
end

