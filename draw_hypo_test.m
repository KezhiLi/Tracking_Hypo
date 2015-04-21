figure, 
%height = Npix_h; 
height = 201; 


iii=ii;
plot(X{iii}.xy(:,1),height-X{iii}.xy(:,2),'b','LineWidth',2);


for jjj=1:15; 
    jjj
    hold on, 
    plot(XX{iii,jjj}.xy(:,1),height-XX{iii,jjj}.xy(:,2),'r'); 
    pause(1)
end; 

plot(XX{iii,49}.xy(:,1),XX{iii,49}.xy(:,2),'g','LineWidth',2);
plot(XX{iii,50}.xy(:,1),XX{iii,50}.xy(:,2),'y','LineWidth',2);

%%%%%%%%%%%%%

hold off


figure, plot(pt_sp_x(:),height-pt_sp_y(:),'b','LineWidth',2)
hold on, plot(X{ii}.xy(:,1),height-X{ii}.xy(:,2),'m','LineWidth',2);


%%%%%%%%%%%

figure, 
height = Npix_h; 

iii=ii;
plot(X{iii}.xy(:,1),X{iii}.xy(:,2),'b','LineWidth',2);


for jjj=1:48; 
    jjj
    hold on, 
    plot(XX{iii,jjj}.xy(:,1),XX{iii,jjj}.xy(:,2),'r'); 
    pause(1)
end; 

%%%%%%%%%%%%%%%%%%%%%%%

figure, 
%height = Npix_h; 
height = 201; 


iii=ii;

for jjj=1:10; 
    jjj
    hold on, 
    plot(X{jjj}.xy(:,1),height-X{jjj}.xy(:,2),'r'); 
    pause(1)
end; 

plot(XX{iii,49}.xy(:,1),XX{iii,49}.xy(:,2),'g','LineWidth',2);
plot(XX{iii,50}.xy(:,1),XX{iii,50}.xy(:,2),'y','LineWidth',2);

%%%%%%%%%%%%%
figure, imshow(C)
%BW_s = edge(II,'sobel'); 
%BW_r = edge(II,'roberts'); 
BW_l = edge(II,'log');   % 
%BW_p = edge(II,'prewitt');
BW_c = edge(II,'canny');  %
%figure,imshow(BW_s)
%figure,imshow(BW_r)
figure,imshow(BW_l)
%figure,imshow(BW_p)
figure,imshow(BW_c)
