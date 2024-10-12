clear;clc;close all;

p = 0.001:0.005:4;
%�����������1000
Maxj = 1000;
%����2��1000 * 1��ȫ0�����������x,y��ֵ���ӿ�����ٶ�
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%����һ��length(a)*2��ȫ������������lamudax��ֵ��lamuday��ֵ���ӿ�����ٶ�
lamuda = zeros(length(p),2);
%����һ��2*1��ȫ������������ÿ��ѭ������lamuda��ֵ���ӿ�����ٶ�
LE  = zeros(2,1);
q=1;

for i = 1:length(p)
    x(1) = 0.6;
    y(1) = 0.7;
    Q = eye(2);
    J0 = [ pi + p(i)*exp(pi)*(pi*y(1) - 1)^2 + 2*q*x(1)*y(1),   q*x(1)^2 + 2*p(i)*pi*exp(pi)*(y(1)*pi - 1)*x(1);
         - q*y(1)^2 + 2*p(i)*pi*exp(pi)*(x(1)*pi - 1)*y(1), pi + p(i)*exp(pi)*(pi*x(1) - 1)^2 - 2*q*x(1)*y(1)];

    B = J0 * Q;
    [Q,R] = qr(B);
    sum = log(abs(diag(R)));
    for j = 1:Maxj-1
        %����new����
        x(j+1) = mod(exp(pi) * (((p(i) * x(j) * (1 - pi*y(j)) ^ 2)))+ q * y(j) * x(j) ^ 2 + pi * x(j),1);
        y(j+1) = mod(exp(pi) * (((p(i) * y(j) * (1 - pi*x(j)) ^ 2)))- q * x(j) * y(j) ^ 2 + pi * y(j),1);
        
        %��ά����ϵͳ��Lyapunovָ��ʹ���ſɱȵ�����ʽ���
        %�ҳ�ƫ������
        J = [ pi + p(i)*exp(pi)*(pi*y(j) - 1)^2 + 2*q*x(j)*y(j),   q*x(j)^2 + 2*p(i)*pi*exp(pi)*(y(j)*pi - 1)*x(j);
         - q*y(j)^2 + 2*p(i)*pi*exp(pi)*(x(j)*pi - 1)*y(j), pi + p(i)*exp(pi)*(pi*x(j) - 1)^2 - 2*q*x(j)*y(j)];
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
plot(p,lamuda(:,1),...
    '-diamond',...
    'color',[144, 44, 252]/255,...
    'MarkerFaceColor',[241, 156, 162]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(p));
hold on ;
plot(p,lamuda(:,2),...
    '-.^',...
    'color',[44, 186, 252]/255,...
    'MarkerFaceColor',[0, 0, 0]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(p));
% y1=line([0,4],[0,0]);
% set(y1,'linestyle','--','color','r','LineWidth',1.4);
lgd=legend('\it{LE_1}','\it{LE_2}','Location','SouthEast'); 
set(lgd,'FontName','Times New Roman','FontSize',17);
set(gca,'LooseInset',get(gca,'TightInset'));
ylabel('LEs','FontName','Times New Roman','FontSize',17);
xlabel('\it{p}','FontName','Times New Roman','FontSize',17);

saveas(gcf,'EPHM_p.png')

