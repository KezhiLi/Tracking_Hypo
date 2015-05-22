xx = texture_63(:,1);
yy = texture_64(:,1);

% close all;
% clear all;
% clc;
% x=input('enter the value of 1st sequence');
% j=input('enter the value of 2nd sequence');

x = xx';
j = yy';
h=fliplr(j);
disp('the 1st sequence is-');
disp(x); 
disp('the 2nd sequence is-');
disp(j);
lx=length(x);
lh=length(h);
n=lx+lh-1;
subplot(3,1,1);
stem(x);
title('1st sequence');
subplot(3,1,2);
stem(j);
title('2nd sequence');
hh=[h zeros(1,n-lh)];
xx=zeros(n);
xx(1:lx,1)=x;
for i=2:n
    for j=2:n
        xx(j,i)=xx(j-1,i-1);
      
    end;
end;
yy=xx*hh';
subplot(3,1,3);
stem(yy);
disp('cross correlate o/p->');
disp(yy');
title('y=cross correlastion of x & j');

