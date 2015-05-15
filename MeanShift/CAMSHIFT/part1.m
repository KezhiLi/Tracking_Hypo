% ----------------------------------------------------------------------
% - CSE 486
% - Project 5
% - Group 8
% - idg101, adi102, jlm522
% ----------------------------------------------------------------------

% Part 1

% cd frames;

for i = 2:10,
    filename1 = sprintf('%3.3i.png', i);
    filename2 = sprintf('%3.3i.png', i-1);
    A = imread(filename1, 'png');
    A = rgb2gray(A);
    
    B = imread(filename2, 'png');
    B = rgb2gray(B);
    
    C = imabsdiff(A,B);
    imwrite(C, sprintf('%3.3i-%3.3i.png', i, i-1), 'png');
end

cd ..

disp('Done.');