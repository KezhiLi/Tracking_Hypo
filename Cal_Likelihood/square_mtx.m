function dis_mtx = square_mtx(row, col, pt, len)
%
% 
% 
% 
% 

A1 = 1:col;
A11 = ones(row,1);
A_c = kron(A11,A1);

A2 = (1: row)';
A22 = ones(1, col);
A_r = kron(A22,A2);

B_r = ones(row,col)*pt(1);
B_c = ones(row,col)*pt(2);

C_1 = abs(A_c - B_r) < len/2;
C_2 = abs(A_r - B_c) < len/2;
dis_mtx = (C_1 + C_2) == 2;


