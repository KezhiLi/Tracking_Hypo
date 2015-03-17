figure, 
height = 175; 

iii=ii;
plot(X{iii}.xy(:,1),X{iii}.xy(:,2),'b','LineWidth',2);


for jjj=1:48; 
    jjj
    hold on, 
    plot(XX{iii,jjj}.xy(:,1),XX{iii,jjj}.xy(:,2),'r'); 
    pause(1)
end; 
plot(XX{iii,49}.xy(:,1),XX{iii,49}.xy(:,2),'g','LineWidth',2);
plot(XX{iii,50}.xy(:,1),XX{iii,50}.xy(:,2),'y','LineWidth',2);



hold off


figure, plot(pt_sp_x(:),height-pt_sp_y(:),'b','LineWidth',2)
hold on, plot(X{ii}.xy(:,1),height-X{ii}.xy(:,2),'m','LineWidth',2);