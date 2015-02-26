function [k, N, Vertices, Lines, Vertices_old] = curvature_areainput(Area)
%%
% Input: BWoutline is a logical matrix decribles the outline of a shape
% Output: k is M x 1 Curvature values
%           N :  Vertices 2 x M
%  Function is written by Kezhi Li of Imperial College London (Nov. 2014)
%%
[mm,nn] = size(Area);
% % record the indexes of points on the boundary
% [fy, fx] = find(BWoutline == 1); %also thresholding of grayscale image can do the same job
% %upside down the image
% fy = mm - fy;

Area = flipud(Area);

% sort the point to make sure that a line links them in a sequence
%data_ordered = image_sort(fx,fy); 
data_ordered = contour_following(Area);
% To have the length of the line sequence
[mm_data_ord, nn_data_ord] = size(data_ordered);
% record the Vertices
Vertices = double([data_ordered(:,2), data_ordered(:,1)]);
% record the Lines link the Vertices. Because they are already sorted, just
% from 1~2,2~3,3~4...
Lines = double([1:mm_data_ord; [2:mm_data_ord 1]]');

% Smooth the Vertices of the line, by applying a local average filter 
Vertices_smooth = zeros(mm_data_ord, nn_data_ord);
% The parameter (0.4, 0.3 etc) can be changed adaptively
for ii=1:mm_data_ord;
    Vertices_smooth(ii,:)= 0.4*Vertices(ii,:)+ 0.3*(Vertices(Lines(find(Lines(:,1)==ii),2),:)+Vertices(find(Lines(:,2)==ii),:));
end
Vertices_old = Vertices;
Vertices = Vertices_smooth;

% Calculate the Curvature
 k=LineCurvature2D(Vertices,Lines);
 figure,  hold on;
 % The curvature direction is calculated.  
 % LineNormal2D_angle3point calculates the curvature based on 2 points that distance 3 pointsaway
 N=LineNormals2D_angle3point(Vertices,Lines);
 % Adjust the magnitude of the green line to make it visible 
 k=k*20;
 % Plot the results
 plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
 %plot([Vertices(:,1) Vertices(:,1)+50.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+50.*N(:,2)]','g');
 plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
 %plot(Vertices(:,1),Vertices(:,2),'r.');
 
 Vertices = round(Vertices);
 
 axis equal;
