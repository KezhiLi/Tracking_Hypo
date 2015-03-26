function dis_mtx = round_mtx(row, col, pt)
%
% 
% 
% 
% 

A1 = 1:col;
A11 = ones(row,1);
A_1 = kron(A11,A1);

A_2 = A_1';

B_1 = ones(row,col)*pt(1);
B_2 = ones(row,col)*pt(2);

dis_mtx = sqrt( (A_1 - B_1).^2 + (A_2 - B_2).^2);
