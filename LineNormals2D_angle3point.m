function N=LineNormals2D_angle3point(Vertices,Lines)
% This function calculates the normals, of the line points
% using the neighbouring points of each contour point, and 
% forward an backward differences on the end points
%
% N=LineNormals2D(V,L)
%
% inputs,
%   V : List of points/vertices 2 x M
% (optional)
%   Lines : A N x 2 list of line pieces, by indices of the vertices
%         (if not set assume Lines=[1 2; 2 3 ; ... ; M-1 M])
%
% outputs,
%   N : The normals of the Vertices 2 x M
%
% Example, Hand
%  load('testdata');
%  N=LineNormals2D(Vertices,Lines);
%  figure,
%  plot([Vertices(:,1) Vertices(:,1)+10*N(:,1)]',[Vertices(:,2) Vertices(:,2)+10*N(:,2)]');
%
% Function is written by modified by Kezhi Li, Nov. 2014
% 

% If no line-indices, assume a x(1) connected with x(2), x(3) with x(4) ...
if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end

% Calculate tangent vectors
DT=Vertices(Lines(:,1),:)-Vertices(Lines(Lines(Lines(Lines(:,2),2),2),2),:);

% Make influence of tangent vector 1/Distance
% (Weighted Central Differences. Points which are closer give a 
% more accurate estimate of the normal)

% LL=sqrt(DT(:,1).^2+DT(:,2).^2);
% DT(:,1)=DT(:,1)./max(LL.^2,eps);
% DT(:,2)=DT(:,2)./max(LL.^2,eps);

D1=zeros(size(Vertices)); D1(Lines(:,1),:)=DT;
D2=zeros(size(Vertices)); D2(Lines(:,2),:)=DT;
D=D1+D2;

mag_D1 = sqrt(D1(:,1).^2+D1(:,2).^2);
mag_D2 = sqrt(D2(:,1).^2+D2(:,2).^2);
% real=aa.*cc-bb.*dd;
% imagi=aa.*dd+bb.*cc;
% mag = sqrt(real.^2+imagi.^2);

mag = D1.*D2;
mag = mag(:,1)+mag(:,2);
theta = acos(mag./(mag_D1.*mag_D2));

% % Normalize the normal
% LL=sqrt(D(:,1).^2+D(:,2).^2);
N(:,1)=-D(:,2);%./LL;
N(:,2)= D(:,1);%./LL;

% keep the direction of angle, and set the magnitude as the angle: theta
mor = sqrt(N(:,1).^2+N(:,2).^2)./theta;
N = N./[mor, mor];
N = [N(end,:);N(1:end-1,:)];
