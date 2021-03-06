function mask_cir = create_cir(C_r, C_c, R, imageSize)
% 
% 
% 
% 
% 
% 
%
R= 4;
L = 2 * R;

% I = ones(L,L);
% imageSize = size(I);
% ci = [R, R, R];     % center and radius of circle ([c_row, c_col, r])
%[xx,yy] = ndgrid((1:imageSize(1))-ci(1)-0.5,(1:imageSize(2))-ci(2)-0.5);

ci = [C_r, C_c, R];
[xx,yy] = ndgrid((1:imageSize(1))-ci(1)-0.5,(1:imageSize(2))-ci(2)-0.5);

mask_cir = (xx.^2 + yy.^2)<ci(3)^2;




