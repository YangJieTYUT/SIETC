clc;clear;close all;

n=20000;

% x(1)=0.6143;y(1)=0.8461;z(1)=0.8461;
% x1(1)=0.6143;y1(1)=0.8461;
% x2(1)=0.6143;y2(1)=0.8461;
% x3(1)=0.6143;y3(1)=0.8461;
% x4(1)=0.6143;y4(1)=0.8461;

x1(1)=0.6;y1(1)=0.7;
x2(1)=0.6;y2(1)=0.7;
x3(1)=0.6;y3(1)=0.7;
x4(1)=0.6;y4(1)=0.7;

for i=1:n-1


% 2D-ILASM
mu = 0.55;rho = 20;
x1(i+1)=sin(pi * mu * (y1(i) + rho - x1(i)) * x1(i) * (1 - x1(i)));
y1(i+1)=sin(pi * mu * (x1(i+1) + rho - y1(i)) * y1(i) * (1 - y1(i)));

% new discrete chaotic map
a = 0.8;b = 1; k = 5;
x2(i+1)= (1-a * x2(i)^2 + y2(i)) * sin(k/x2(i));
y2(i+1)= b* x2(i)* sin(k/x2(i));

% 2D-SLTCM
alpha = 1;w=0.2;
if y3(i)< 0.5
    x3(i+1)=alpha * sin(pi * ((3 + (4 - w) * y3(i) * 0.5) * x3(i) * (1 - x3(i))));
else
    x3(i+1)=alpha * sin(pi * ((3 + (4 - w) * (1 - y3(i)) * 0.5) * x3(i) * (1 - x3(i))));
end
if x3(i+1)< 0.5
    y3(i+1)=alpha * sin(pi * ((3 + (4 - w) * x3(i+1) * 0.5) * y3(i) * (1 - y3(i))));
else
    y3(i+1)=alpha * sin(pi * ((3 + (4 - w) * (1 - x3(i+1)) * 0.5) * y3(i) * (1 - y3(i))));
end


% 2D-EPHM
p = 1;q = 1;
x4(i + 1) = mod(exp(pi) * (((p * x4(i) * (1 - pi*y4(i)) ^ 2)))+ q * y4(i) * x4(i) ^ 2 + pi * x4(i),1);
y4(i + 1) = mod(exp(pi) * (((p * y4(i) * (1 - pi*x4(i)) ^ 2)))- q * x4(i) * y4(i) ^ 2 + pi * y4(i),1);
end
%% 2D»­Í¼
figure;
plot(x2,y2,'b.','markersize',0.1,'HandleVisibility', 'off');hold on;
scatter(x1(1),y1(1),50,'r','filled','pentagram');
lgd=legend('(\itx\rm_{0},\ity\rm_{0})','Location','southeast');
set(lgd,'FontName','Times New Roman','FontSize',15);
% xticks([-1,-0.5,0,0.5,1])
% yticks([-1,-0.5,0,0.5,1])
ylabel('\it{y_i}','FontName','Times New Roman','FontSize',15);
xlabel('\it{x_i}','FontName','Times New Roman','FontSize',15);
set(gca,'FontSize',15,'FontName','Times New Roman');
set(gca,'LooseInset',get(gca,'TightInset'));
