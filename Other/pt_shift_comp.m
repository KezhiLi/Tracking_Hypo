function X = pt_shift_comp(X, CMs, varargin)
% X = pt_shift_comp(X, CMs, k)
% function to adjust coordinates of skeleton points according to the point
% shift record matrix CMs
% 
% Input: X:  a cell that records all skeleton points information
%        CMs: an 2*M matrix that records the coordinates of the central
%        point of the cropped video. M is the number of frames. If the data of CMs in one column    
%        changes to the next column, it means that the 'cemera' shifts
%        k: the index of frame
%        k_sure: the index of last good frame
% Output: X, with new X{ii},xy
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 12/06/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.


% for all X{ii}, ii = 1, ...
for ii = 1:size(X,2);
    % find the central point change, switch x,y coordinates
    if nargin == 3
        k = varargin{1};
        shift1 = fliplr((double(CMs(:,k) - CMs(:,k-1)))');
    elseif nargin == 4
        k = varargin{1};  
        k_sure = varargin{2};
        shift1 = fliplr((double(CMs(:,k) - CMs(:,k_sure)))');
    end
    % if the shift is not zero
    if sum(abs(shift1)) > 0;
        % create a tall matrix that corresponds for x,y coordinates,
        % respectively 
        CM_shift = kron(shift1,ones(size(X{ii}.xy,1),1));
        % adjust the coordinates of X{ii}.xy
        X{ii}.xy = X{ii}.xy - CM_shift;
    end
    
end