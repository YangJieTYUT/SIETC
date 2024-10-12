clear;clc;close all;

b = 0.001:0.005:4;
%�����������1000
Maxj = 1000;
%����2��1000 * 1��ȫ0�����������x,y��ֵ���ӿ�����ٶ�
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%����һ��length(a)*2��ȫ������������lamudax��ֵ��lamuday��ֵ���ӿ�����ٶ�
lamuda = zeros(length(b),2);
%����һ��2*1��ȫ������������ÿ��ѭ������lamuda��ֵ���ӿ�����ٶ�
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
        %����new����

        x(j+1)= (1-a * x(j)^2 + y(j)) * sin(k/x(j));
        y(j+1)= b(i)* x(j)* sin(k/x(j));
        
        %��ά����ϵͳ��Lyapunovָ��ʹ���ſɱȵ�����ʽ���
        %�ҳ�ƫ������
        J = [ - 2*a*x(j)*sin(k/x(j)) - (k*cos(k/(1))*(- a*x(j)^2 + y(j) + 1))/(x(j)^2), sin(k/x(j));
        b(i)*sin(k/x(j)) - (b(i)*k*cos(k/x(j)))/x(j),        0];
        %�˵�λ�����ԭ�������Ϊ�˼ӿ�����ٶȣ���ʼ����һ��2*2�ľ���
        B = J * Q;
        %����qr�ֽ⣬��B����ֽ����������Q�������Ǿ���R��R��������Խ��ߵ�ֵ��Ϊlamuda��ֵ
        [Q,R] = qr(B);
        %��ά����ϵͳ��Lyapunovָ��1/n*sum(ln(lumada));
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

