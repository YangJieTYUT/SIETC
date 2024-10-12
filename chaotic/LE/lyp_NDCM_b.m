clear;clc;close all;

b = 0.001:0.005:4;
%定义迭代次数1000
Maxj = 1000;
%定义2个1000 * 1的全0矩阵用来存放x,y的值，加快迭代速度
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%定义一个length(a)*2的全零矩阵用来存放lamudax的值、lamuday的值，加快迭代速度
lamuda = zeros(length(b),2);
%定义一个2*1的全零矩阵用来存放每次循环两个lamuda的值，加快迭代速度
LE  = zeros(2,1);
        a = 0.8;
%         b = 1; 
        k = 0.5;

for i = 1:length(b)
    x(1) = 0.6;
    y(1) = 0.7;
    Q = eye(2);
    J0 = [ - 2*a*x(1)*sin(k/x(1)) - (k*cos(k/(1))*(- a*x(1)^2 + y(1) + 1))/(x(1)^2), sin(k/x(1));
        b(i)*sin(k/x(1)) - (b(i)*k*cos(k/x(1)))/x(1),        0];

    B = J0 * Q;
    [Q,R] = qr(B);
    sum = log(abs(diag(R)));
    for j = 1:Maxj-1
        %定义new混沌

        x(j+1)= (1-a * x(j)^2 + y(j)) * sin(k/x(j));
        y(j+1)= b(i)* x(j)* sin(k/x(j));
        
        %高维混沌系统的Lyapunov指数使用雅可比迭代方式求解
        %找出偏导矩阵
        J = [ - 2*a*x(j)*sin(k/x(j)) - (k*cos(k/(1))*(- a*x(j)^2 + y(j) + 1))/(x(j)^2), sin(k/x(j));
        b(i)*sin(k/x(j)) - (b(i)*k*cos(k/x(j)))/x(j),        0];
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
plot(b,lamuda(:,1),...
    '-diamond',...
    'color',[144, 44, 252]/255,...
    'MarkerFaceColor',[241, 156, 162]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(b));
hold on ;
plot(b,lamuda(:,2),...
    '-.^',...
    'color',[44, 186, 252]/255,...
    'MarkerFaceColor',[0, 0, 0]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(b));
y1=line([0,4],[0,0]);
set(y1,'linestyle','--','color','r','LineWidth',1.4);
lgd=legend('\it{LE_1}','\it{LE_2}','Location','SouthEast'); 
set(lgd,'FontName','Times New Roman','FontSize',17);
set(gca,'LooseInset',get(gca,'TightInset'));
ylabel('LEs','FontName','Times New Roman','FontSize',17);
xlabel('\it{b}','FontName','Times New Roman','FontSize',17);

saveas(gcf,'NDCM_b.png')

