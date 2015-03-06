function fxys = staightData(img, polyDegree)
%make an input membrane black and white image straight
%Syntax: xy = staightData(img, polyDegree)
%outptu xy is (N, 2) matrix
%polyDegree is optional degree of fitted polynom, default 5

PLOT_OUTPUTS = 1;
%note that plotted fy values are up-and-down flipped with comparison to the input image

if nargin < 2
    polyDegree = 5;
end

%step 1 - get x and y coordinates
[fy, fx] = find(img == 1); %also thresholding of grayscale image can do the same job
if PLOT_OUTPUTS
   figure, plot(fx, fy, '.') 
   title('Input pixels')
end

%step 2 - rotate image to be "left to right" horizontal - to ensure successful polyfit
c = cov([fx fy]);
[V, L] = eig(c);
[~, idx] = sort(diag(L), 'descend');
fxy = [fx fy] * V(:, idx);
if PLOT_OUTPUTS
   figure, plot(fxy(:, 1), fxy(:,2), 'k.') 
   title('Horizontal rotation of input data')
end

%step 3 - polynomial fit of the data
p = polyfit(fxy(:,1), fxy(:,2), polyDegree);
if PLOT_OUTPUTS
    xx = min(fxy(:,1)):max(fxy(:,1));
    yy = polyval(p, xx);
    figure, plot(fxy(:, 1), fxy(:,2), '.')
    title(sprintf('Polynomial fit (degree %d)', polyDegree));
    hold on
    plot(xx, yy, 'r')
    hold off
end

%step 4 - subtration of the polyfit
fxys = fxy;
fxys(:, 2) = fxys(:, 2) - polyval(p, fxy(:, 1));
if PLOT_OUTPUTS
    figure, plot(fxys(:, 1), fxys(:,2), 'r.')
    title('Final wiggles after subtraction of the polynomial fit')
end
