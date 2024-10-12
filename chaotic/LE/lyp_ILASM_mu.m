clear;clc;close all;

mu = 0.001:0.005:4;
%定义迭代次数1000
Maxj = 1000;
%定义2个1000 * 1的全0矩阵用来存放x,y的值，加快迭代速度
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%定义一个length(a)*2的全零矩阵用来存放lamudax的值、lamuday的值，加快迭代速度
lamuda = zeros(length(mu),2);
%定义一个2*1的全零矩阵用来存放每次循环两个lamuda的值，加快迭代速度
LE  = zeros(2,1);
rho  = 4;

for i = 1:length(mu)
    x(1) = 0.2;
    y(1) = 0.3;
    Q = eye(2);
    J0 = [ -cos(mu(i)*x(1)*pi*(x(1) - 1)*(rho - x(1) + y(1)))*(mu(i)*pi*(x(1) - 1)*(rho - x(1) + y(1)) - pi*mu(i)*x(1)*(x(1) - 1) + mu(i)*x(1)*pi*(rho - x(1) + y(1))), -mu(i)*x(1)*pi*cos(mu(i)*x(1)*pi*(x(1) - 1)*(rho - x(1) + y(1)))*(x(1) - 1);
       -mu(i)*y(1)*pi*cos(mu(i)*y(1)*pi*(y(1) - 1)*(rho + x(1) - y(1)))*(y(1) - 1), -cos(mu(i)*y(1)*pi*(y(1) - 1)*(rho + x(1) - y(1)))*(mu(i)*pi*(y(1) - 1)*(rho + x(1) - y(1)) - pi*mu(i)*y(1)*(y(1) - 1) + mu(i)*y(1)*pi*(rho + x(1) - y(1)))];
    B = J0 * Q;
    [Q,R] = qr(B);
    sum = log(abs(diag(R)));
    for j = 1:Maxj-1
        %定义new混沌
        x(j+1)=sin(pi * mu(i) * (y(j) + rho - x(j)) * x(j) * (1 - x(j)));
        y(j+1)=sin(pi * mu(i) * (x(j+1) + rho - y(j)) * y(j) * (1 - y(j)));
        
        %高维混沌系统的Lyapunov指数使用雅可比迭代方式求解
        %找出偏导矩阵
        J = [ -cos(mu(i)*x(j)*pi*(x(j) - 1)*(rho - x(j) + y(j)))*(mu(i)*pi*(x(j) - 1)*(rho - x(j) + y(j)) - pi*mu(i)*x(j)*(x(j) - 1) + mu(i)*x(j)*pi*(rho - x(j) + y(j))), -mu(i)*x(j)*pi*cos(mu(i)*x(j)*pi*(x(j) - 1)*(rho - x(j) + y(j)))*(x(j) - 1);
              -mu(i)*y(j)*pi*cos(mu(i)*y(j)*pi*(y(j) - 1)*(rho + x(j) - y(j)))*(y(j) - 1), -cos(mu(i)*y(j)*pi*(y(j) - 1)*(rho + x(j) - y(j)))*(mu(i)*pi*(y(j) - 1)*(rho + x(j) - y(j)) - pi*mu(i)*y(j)*(y(j) - 1) + mu(i)*y(j)*pi*(rho + x(j) - y(j)))];
        %乘单位矩阵的原因可能是为了加快迭代速度？初始化了一个2*2的矩阵
        B = J * Q;
        %进行qr分解，将B矩阵分解成正交矩阵Q和上三角矩阵R，R矩阵的正对角线的值即为lamuda的值
        [Q,R] = qr(B);
        %高维混沌系统的Lyapunov指数1/n*sum(ln(lumada));
        sum = sum + log(abs(diag(R)));
        
    end
    LE = sum / Maxj;
    lamuda(i,1) = LE(1);
    lamuda(i,2) = LE(2);
    
end
% plot(mu,lamuda(:,1),'g-','LineWidth',1);
% hold on ;
% plot(mu,lamuda(:,2),'b-','LineWidth',1);
%%
figure;
plot(mu,lamuda(:,1),...
    '-diamond',...
    'color',[144, 44, 252]/255,...
    'MarkerFaceColor',[241, 156, 162]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(mu));
hold on ;
plot(mu,lamuda(:,2),...
    '-.^',...
    'color',[44, 186, 252]/255,...
    'MarkerFaceColor',[0, 0, 0]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(mu));
y1=line([0,4],[0,0]);
set(y1,'linestyle','--','color','r','LineWidth',1.4);
lgd=legend('\it{LE_1}','\it{LE_2}','Location','SouthEast'); 
set(lgd,'FontName','Times New Roman','FontSize',17);
set(gca,'LooseInset',get(gca,'TightInset'));
ylabel('LEs','FontName','Times New Roman','FontSize',17);
xlabel('\mu','FontName','Times New Roman','FontSize',17);

saveas(gcf,'ILASM_mu.png')

