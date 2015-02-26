function show_better_estimated_Worm_body_hypo(worm_show, Y_k, width)

%X_mean = mean(X, 2);

figure(1)
image(Y_k)
title('+++ Showing nearest point +++')

    seg_len = 8;
    m_fre_pt = size(worm_show,1);
    t = 1:m_fre_pt ;
    ts = 1:1/(seg_len*2):m_fre_pt;
    
    tc = 1:2*m_fre_pt ;
    tsc = 1:1/(seg_len*2):2*m_fre_pt;


worm_ske1 = worm_show;
worm_ske_curv_x =  0.7*interp1(t,worm_ske1(:,1),ts,'spline')+0.3*interp1(t,worm_ske1(:,1),ts,'linear') ;
worm_ske_curv_y =  0.7*interp1(t,worm_ske1(:,2),ts,'spline')+0.3*interp1(t,worm_ske1(:,2),ts,'linear') ;

[worm_ske1_T,worm_ske1_N] = frenet_TN(worm_ske1(:,1),worm_ske1(:,2));
[worm_shape,worm_body] = ske2shape(worm_ske1, worm_ske1_N, width, -0.4);
worm_shape1 = round(worm_shape);
worm_shape_x =  round(0.7*interp1(tc,worm_shape1(:,1),tsc,'spline')+0.3*interp1(tc,worm_shape1(:,1),tsc,'linear') );
worm_shape_y =  round(0.7*interp1(tc,worm_shape1(:,2),tsc,'spline')+0.3*interp1(tc,worm_shape1(:,2),tsc,'linear') );


% worm_shape =ske2shape(X{1}.xy, X{1}.N, width, -0.25);
% tt=5;
% for ii=2:tt;
% worm_shape = worm_shape+ske2shape(X{ii}.xy, X{1}.N, width, -0.25);
% end
% worm_shape = round(worm_shape/tt);
% 
%     seg_len = 8;
%     m_fre_pt = size(X{1,1}.xy,1);
%     t = 1:2*m_fre_pt ;
%     ts = 1:1/(seg_len*2):2*m_fre_pt;
% worm_shape1 = round(worm_shape);
% worm_shape_x =  round(0.5*interp1(t,worm_shape1(:,1),ts,'spline')+0.5*interp1(t,worm_shape1(:,1),ts,'linear') );
% worm_shape_y =  round(0.5*interp1(t,worm_shape1(:,2),ts,'spline')+0.5*interp1(t,worm_shape1(:,2),ts,'linear') );


hold on


 plot(worm_shape_x,worm_shape_y,'LineWidth',2,'Color',[1 0.1 0.1]);
plot(worm_ske_curv_x,worm_ske_curv_y,'LineWidth',2,'Color',[1 0.1 0.1]);

% plot(X{1}.xy(:,1), X{1}.xy(:,2) ,'b*');
% plot(worm_shape1(:,1),worm_shape1(:,2),'LineWidth',2,'Color',[1 0.1 0.1]);
% 
% plot(X{2}.xy(:,1), X{2}.xy(:,2) ,'b*');
% plot(worm_shape2(:,1),worm_shape2(:,2),'LineWidth',2,'Color',[1 0.1 0.1]);
% 
% plot(X{3}.xy(:,1), X{3}.xy(:,2) ,'b*');
% plot(worm_shape3(:,1),worm_shape3(:,2),'LineWidth',2,'Color',[1 0.1 0.1]);

hold off 


% num_shown = 50;
% 
% index = find(Q==max(Q));
% 
% 
% hold on
% 
% %plot(X_mean(2,:), X_mean(1,:), 'h', 'MarkerSize', 16, 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y')
% plot(mean(X(1,index)), mean(X(2,index)), 'h', 'MarkerSize', 16, 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y')
% 
% plot(X(1,1:num_shown),X(2,1:num_shown),'b*')
% hold off









%drawnow
