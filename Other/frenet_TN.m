function [TT,NN] = frenet_TN(x,y,z)
% FRENET_TN - Frenet-Serret Space Curve Invarients
% A modified simpler version of Frenet, to accelerate the process
%   
%   [T,N] = frenet_TN(x,y);
%   [T,N] = frenet_TN(x,y,z);
% 
%   Returns the 3 vector and 2 scaler invarients of a space curve defined
%   by vectors x,y and z.  If z is omitted then the curve is only a 2D,
%   but the equations are still valid.
% 
%    _    r'
%    T = ----  (Tangent)
%        |r'|
% 
%    _    T'
%    N = ----  (Normal)
%        |T'|
%    _   _   _
%    B = T x N (Binormal)
% 
%    k = |T'|  (Curvature)
% 
%    t = dot(-B',N) (Torsion)
% 
% 
%    Example:
%    theta = 2*pi*linspace(0,2,100);
%    x = cos(theta);
%    y = sin(theta);
%    z = theta/(2*pi);
%    [T,N,B,k,t] = frenet(x,y,z);
%    line(x,y,z), hold on
%    quiver3(x,y,z,T(:,1),T(:,2),T(:,3),'color','r')
%    quiver3(x,y,z,N(:,1),N(:,2),N(:,3),'color','g')
%    quiver3(x,y,z,B(:,1),B(:,2),B(:,3),'color','b')
%    legend('Curve','Tangent','Normal','Binormal')
% 
% 
% See also: GRADIENT
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% modified function frenet.m to make it faster
% 18/01/2015


if nargin == 2,
    z = zeros(size(x));
end

% CONVERT TO COLUMN VECTOR
x = x(:);
y = y(:);
z = z(:);

% SPEED OF CURVE
% dx = gradient(x);
% dy = gradient(y);
% dz = gradient(z);
% dr = [dx dy dz];

dx1 = DGradient(x,length(x),1, '1stOrder');
dy1 = DGradient(y,length(y),1, '1stOrder');
dz1 = DGradient(z,length(z),1, '1stOrder');
dr1 = [dx1 dy1 dz1];


% TANGENT
% TT = dr;
% T = dr./mag(dr,3);

TT = dr1;
%T1 = dr1./mag(dr1,3);

% DERIVIATIVE OF TANGENT
% dTx =  gradient(T(:,1));
% dTy =  gradient(T(:,2));
% dTz =  gradient(T(:,3));
% dT = [dTx dTy dTz];

size_T1 = size(TT,1);
dTx1 = DGradient(TT(:,1),size_T1,1, '1stOrder');
dTy1 = DGradient(TT(:,2),size_T1,1, '1stOrder');
dTz1 = DGradient(TT(:,3),size_T1,1, '1stOrder');
dT1 = [dTx1 dTy1 dTz1];
dT11 =  [dTx1 dTy1];

% NORMAL
NN = dT11./mag(dT1,nargin);


% N = dT./mag(dT,3);
% % BINORMAL
% B = cross(T,N);
% BB = cross(TT,NN);
% % CURVATURE
% % k = mag(dT,1);
% k = mag(cross(dr,ddr),1)./((mag(dr,1)).^3);
% % TORSION
% t = dot(-B,N,2);




function N = mag(T,n)
% MAGNATUDE OF A VECTOR (Nx3)
%  M = mag(U)
N = sum(abs(T).^2,2).^(1/2);
d = find(N==0); 
N(d) = eps*ones(size(d));
N = N(:,ones(n,1));