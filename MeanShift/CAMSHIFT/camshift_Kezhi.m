% ----------------------------------------------------------------------
% - CSE 486
% - Project 5
% - Group 8
% - idg101, adi102, jlm522
% ----------------------------------------------------------------------

clc;
disp('Running...');
close all;
clear;
%cd frames;

addpath(genpath('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_16\MeanShift\CAMSHIFT\.'));


% ----------------------------------------------------------------------
% Constants
    WHITE = 255;
 
    % Threshold of convergence (in pixels for each deg. of freedom)
	T = 1;
    
    % Number of pixels to expand search window by.
    P = 5;

% ----------------------------------------------------------------------
% Mean Shift

avi = VideoWriter('output1.avi');
%avi = avifile('output.avi');

open(avi);

% Initial search window size.
W = [10 10];
%W = [80 94];
%W = [200 350];

% Initial location of search window.
%L = [95 193];
L = [102, 121];  % [91, 81]   [ 106 85]

% For plotting motion
Xpoints=[];
Ypoints=[];

disp('Frame: Coordinates');

for frame = 40:85,
	filename = sprintf('%3.3i.jpg', frame);
	R = imread(filename, 'jpg');
	R = 255- R;
	% Convert image from RGB space to HSV space
    % R_rgb = gray2rgb(R);
    % I_hsv = rgb2hsv(R_rgb);
	
	% Extract the hue information
 	% I1 = I_hsv(:,:,3);
    % I = roicolor(I1, 0.83, 1.0);
    %I = double(R);
     I = roicolor(R, 70,255);
    
    % Initialization
    oldCamL = [0 0];
    MeanConverging = 1;
% ----------------------------------------------------------------------
    % Create search window box on image.
    for i = L(1) : L(1)+W(1),
        x = i;
        y = L(2);
        if x > size(I,1) | y > size(I,2) | x < 1 | y < 1
            continue;
        else
            R(x, y,:) = 0;
        end
    end
    
    for i = L(1) : L(1)+W(1),
        x = i;
        y = L(2) + W(2);
        if x > size(I,1) | y > size(I,2) | x < 1 | y < 1
            continue;
        else
            R(x, y, :) = 0;
        end
    end    
    
    for i = L(2) : L(2)+W(2),
        x = L(1);
        y = i;
        if x > size(I,1) | y > size(I,2) | x < 1 | y < 1
            continue;
        else
            R(x, y, :) = 0;
        end
    end    

    for i = L(2) : L(2)+W(2),
        x = L(1)+W(1);
        y = i;
        if x > size(I,1) | y > size(I,2) | x < 1 | y < 1
            continue;
        else
            R(x, y, :) = 0;
        end
    end    
% ----------------------------------------------------------------------  
    tt = 0;
    while MeanConverging && tt < 1000,
        tt = tt + 1

        % Compute centroid of search window
		M00 = 0.0;
		for i = L(1)-P : (L(1)+W(1)+P),
            for j = L(2)-P : (L(2)+W(2)+P),
                if i > size(I,1) | j > size(I,2) | i < 1 | j < 1
                    continue;
                end
                M00 = M00 + double(I(i,j));
            end
		end
		
		M10 = 0.0;
		for i = L(1)-P : (L(1)+W(1)+P),
            for j = L(2)-P : (L(2)+W(2)+P),
                if i > size(I,1) | j > size(I,2) | i < 1 | j < 1
                    continue;
                end
                M10 = M10 + i * double(I(i,j));
            end
		end
		
		M01 = 0.0;
		for i = L(1)-P : (L(1)+W(1)+P),
            for j = L(2)-P : (L(2)+W(2)+P),
                if i > size(I,1) | j > size(I,2)| i < 1 | j < 1
                    continue;
                end                
                M01 = M01 + j * double(I(i,j));
            end
        end
        
        if M10 == 0 && M00 == 0;
            xc = L(1);
            yc = L(2); 
        else
            xc = round(M10 / M00);
            yc = round(M01 / M00);
        end

		oldL = L;
		L = [floor(xc - (W(1)/2)) floor(yc - (W(2)/2))];
       
        % Check threshold
        if abs(oldL(1)-L(1)) < T && abs(oldL(2)-L(2)) < T
            MeanConverging = 0;
        end
    end
    
    % We now know the centroid and M00.
    % This information is used to alter the search window size.
    
%     % Adjust window size
%     % s = round(1.1 * sqrt(M00));
%     s = round(1.1 * sqrt(M00/255));
%     W = [ s      floor(1.2*s) ];
      L = [floor(xc - (W(1)/2)) floor(yc - (W(2)/2))];
    
    % Output the centroid's coordinates
    disp(sprintf('%3i:   %3i, %3i', frame, xc, yc));
    Xpoints = [Xpoints xc];
    Ypoints = [Ypoints yc];
    
    % Superimpose plus sign on to centroid of hand.
    plus_sign_mask = [0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     1 1 1 1 1 1 1 1 1 1 1 1 1;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0 0 0 0 0 0];
	sizeM = size(plus_sign_mask);
	for i = -floor(sizeM(1) / 2):floor(sizeM(1) / 2),
        for j = -floor(sizeM(2) / 2):floor(sizeM(2) / 2),
            if plus_sign_mask(i+1+floor(sizeM(1) / 2), j+1+floor(sizeM(2) / 2)) == 1
                R(i+xc, j+yc, :) = WHITE;
            end
        end
	end	
% ----------------------------------------------------------------------
    % Display the probability image.
    %I = rgb2hsv(R_rgb);
    I = double(R);
    imshow(R)
    
%     S = [];
%     S(:,:,1) = I(:,:,1);
%     S(:,:,2) = I(:,:,1);    
%     S(:,:,3) = I(:,:,1);
    
	% Extract the hue information
    % avi = addframe(avi, S);
    fra = getframe;
    writeVideo(avi, fra);
    pause(0.5)
end

disp('AVI move parameters:');
close(avi);
% ----------------------------------------------------------------------
figure, plot(Ypoints,Xpoints, 'go' , Ypoints, Xpoints);
axis([0 320 0 240]);


% ----------------------------------------------------------------------
cd ..
disp('Done.');
% ----------------------------------------------------------------------
