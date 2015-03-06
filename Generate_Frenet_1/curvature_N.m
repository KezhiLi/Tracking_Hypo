function [N, Vertices, Lines, Vertices_orignal] = curvature_N(BWoutline)
%%
% Input: BWoutline is a logical matrix decribles the outline of a shape
% Output:   N :  Vertors M * 2 matrix describle the curvature
%           Vertices: M *2 matrix after smoothing desribles the Vertices
%           Lines: M *2 matrix desribles the lines between Vertices
%           Vertices_original: original M*2 matrix before smoothing 
%  Function is written by Kezhi Li of Imperial College London (Nov. 2014)
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%%
[mm,nn] = size(BWoutline);
% record the indexes of points on the boundary
[fy, fx] = find(BWoutline == 1); %also thresholding of grayscale image can do the same job
%upside down the image
fy = mm - fy;
% sort the point to make sure that a line links them in a sequence
data_ordered = image_sort(fx,fy); 
% To have the length of the line sequence
[mm_data_ord, nn_data_ord] = size(data_ordered);
% record the Vertices
Vertices = double(data_ordered);
% record the Lines link the Vertices. Because they are already sorted, just
% from 1~2,2~3,3~4...
Lines = double([1:mm_data_ord; [2:mm_data_ord 1]]');

% Smooth the Vertices of the line, by applying a local average filter 
Vertices_smooth = zeros(mm_data_ord, nn_data_ord);
% The parameter (0.4, 0.3 etc) can be changed adaptively
for ii=1:mm_data_ord;
    Vertices_smooth(ii,:)= 0.4*Vertices(ii,:)+ 0.3*(Vertices(Lines(find(Lines(:,1)==ii),2),:)+Vertices(find(Lines(:,2)==ii),:));
end
Vertices_orignal = Vertices;
Vertices = Vertices_smooth;

% Calculate the Curvature
% k=LineCurvature2D(Vertices,Lines);

 % The curvature direction is calculated.  
 % LineNormal2D_angle3point calculates the curvature based on 2 points that distance 3 pointsaway
 N=LineNormals2D_angle3point(Vertices,Lines);

