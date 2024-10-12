clear;clc;close all;

mu = 0.001:0.005:4;
%�����������1000
Maxj = 1000;
%����2��1000 * 1��ȫ0�����������x,y��ֵ���ӿ�����ٶ�
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%����һ��length(a)*2��ȫ������������lamudax��ֵ��lamuday��ֵ���ӿ�����ٶ�
lamuda = zeros(length(mu),2);
%����һ��2*1��ȫ������������ÿ��ѭ������lamuda��ֵ���ӿ�����ٶ�
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
        %����new����
        x(j+1)=sin(pi * mu(i) * (y(j) + rho - x(j)) * x(j) * (1 - x(j)));
        y(j+1)=sin(pi * mu(i) * (x(j+1) + rho - y(j)) * y(j) * (1 - y(j)));
        
        %��ά����ϵͳ��Lyapunovָ��ʹ���ſɱȵ�����ʽ���
        %�ҳ�ƫ������
        J = [ -cos(mu(i)*x(j)*pi*(x(j) - 1)*(rho - x(j) + y(j)))*(mu(i)*pi*(x(j) - 1)*(rho - x(j) + y(j)) - pi*mu(i)*x(j)*(x(j) - 1) + mu(i)*x(j)*pi*(rho - x(j) + y(j))), -mu(i)*x(j)*pi*cos(mu(i)*x(j)*pi*(x(j) - 1)*(rho - x(j) + y(j)))*(x(j) - 1);
              -mu(i)*y(j)*pi*cos(mu(i)*y(j)*pi*(y(j) - 1)*(rho + x(j) - y(j)))*(y(j) - 1), -cos(mu(i)*y(j)*pi*(y(j) - 1)*(rho + x(j) - y(j)))*(mu(i)*pi*(y(j) - 1)*(rho + x(j) - y(j)) - pi*mu(i)*y(j)*(y(j) - 1) + mu(i)*y(j)*pi*(rho + x(j) - y(j)))];
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

