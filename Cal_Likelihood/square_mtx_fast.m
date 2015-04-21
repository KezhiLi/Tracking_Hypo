function dis_mtx = square_mtx_fast(dis_mtx, pt, len)
% add a square matrix mask in the matrix 'dis_mtx'
% 
% Input: dis_mtx: an M*N matrix, the input matrix
%        pt: a 1*2 vector, central point location of the sqaure
%        len: the side length of the square matrix
% Output: dis_mtx: an M*N matrix, with a squared added in
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 16/04/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


Y_1 = size(dis_mtx, 1);
Y_2 = size(dis_mtx, 2);


left = round(max(pt(1)-len,1));
right =round(min(pt(1)+len,Y_2));
bottom = round(max(pt(2)-len,1));
top = round(min(pt(2)+len,Y_1));

dis_mtx(bottom:top,left:right) = 1.5;


