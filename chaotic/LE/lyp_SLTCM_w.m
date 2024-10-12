clear;clc;close all;

w = 0.001:0.01:4;
%�����������1000
Maxj = 1000;
%����2��1000 * 1��ȫ0�����������x,y��ֵ���ӿ�����ٶ�
x = zeros(Maxj,1);
y = zeros(Maxj,1);
%����һ��length(a)*2��ȫ������������lamudax��ֵ��lamuday��ֵ���ӿ�����ٶ�
lamuda = zeros(length(w),2);
%����һ��2*1��ȫ������������ÿ��ѭ������lamuda��ֵ���ӿ�����ٶ�
LE  = zeros(2,1);

alpha = 1;
for i = 1:length(w)
    x(1) = 0.6;
    y(1) = 0.7;
    Q = eye(2);
        J0 =[ -alpha*cos(x(1)*pi*(((w(i) - 4)*(y(1) - 1))/2 + 3)*(x(1) - 1))*(x(1)*pi*(((w(i) - 4)*(y(1) - 1))/2 + 3) + pi*(((w(i) - 4)*(y(1) - 1))/2 + 3)*(x(1) - 1)),   -alpha*x(1)*pi*cos(x(1)*pi*(((w(i) - 4)*(y(1) - 1))/2 + 3)*(x(1) - 1))*(w(i)/2 - 2)*(x(1) - 1);
            -alpha*y(1)*pi*cos(y(1)*pi*(((w(i) - 4)*(x(1) - 1))/2 + 3)*(y(1) - 1))*(w(i)/2 - 2)*(y(1) - 1), -alpha*cos(y(1)*pi*(((w(i) - 4)*(x(1) - 1))/2 + 3)*(y(1) - 1))*(y(1)*pi*(((w(i) - 4)*(x(1) - 1))/2 + 3) + pi*(((w(i) - 4)*(x(1) - 1))/2 + 3)*(y(1) - 1))];
    
    
        B = J0 * Q;
        [Q,R] = qr(B);
        sum = log(abs(diag(R)));
%     sum=0;
    for j = 1:Maxj-1
        %����new����
        if y(j)< 0.5
            x(j+1)=alpha * sin(pi * ((3 + (4 - w(i)) * y(j) * 0.5) * x(j) * (1 - x(j))));
        else
            x(j+1)=alpha * sin(pi * ((3 + (4 - w(i)) * (1 - y(j)) * 0.5) * x(j) * (1 - x(j))));
        end
        if x(j+1)< 0.5
            y(j+1)=alpha * sin(pi * ((3 + (4 - w(i)) * x(j+1) * 0.5) * y(j) * (1 - y(j))));
        else
            y(j+1)=alpha * sin(pi * ((3 + (4 - w(i)) * (1 - x(j+1)) * 0.5) * y(j) * (1 - y(j))));
        end
        %��ά����ϵͳ��Lyapunovָ��ʹ���ſɱȵ�����ʽ���
        %�ҳ�ƫ������
%         if y(j)< 0.5 && x(j+1)< 0.5
%             J=[ alpha*cos(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1))*(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3) + pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1)), alpha*x(j+1)*pi*cos(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1))*(w(i)/2 - 2)*(x(j+1) - 1);
%                 alpha*y(j)*pi*cos(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), alpha*cos(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3) + pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))];
%         elseif y(j)< 0.5 && x(j+1)>= 0.5
%             J=[ alpha*cos(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1))*(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3) + pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1)), alpha*x(j+1)*pi*cos(x(j+1)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j+1) - 1))*(w(i)/2 - 2)*(x(j+1) - 1);
%                 -alpha*y(j)*pi*cos(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), -alpha*cos(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))*(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3) + pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))];
%         elseif y(j)>= 0.5 && x(j+1)< 0.5
%             J=[ -alpha*cos(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1))*(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1)),   -alpha*x(j+1)*pi*cos(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1))*(w(i)/2 - 2)*(x(j+1) - 1);
%                 alpha*y(j)*pi*cos(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), alpha*cos(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(y(j)*pi*((x(j+1)*(w(i) - 4))/2 - 3) + pi*((x(j+1)*(w(i) - 4))/2 - 3)*(y(j) - 1))];
%         elseif y(j)>= 0.5 && x(j+1)>= 0.5
%             J=[ -alpha*cos(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1))*(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1)),   -alpha*x(j+1)*pi*cos(x(j+1)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j+1) - 1))*(w(i)/2 - 2)*(x(j+1) - 1);
%                 -alpha*y(j)*pi*cos(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), -alpha*cos(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))*(y(j)*pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3) + pi*(((w(i) - 4)*(x(j+1) - 1))/2 + 3)*(y(j) - 1))];
%         end
        %j
        if y(j)< 0.5 && x(j+1)< 0.5
            J=[ alpha*cos(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1))*(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3) + pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1)), alpha*x(j)*pi*cos(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1))*(w(i)/2 - 2)*(x(j) - 1);
                alpha*y(j)*pi*cos(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), alpha*cos(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3) + pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))];
        elseif y(j)< 0.5 && x(j+1)>= 0.5
            J=[ alpha*cos(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1))*(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3) + pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1)), alpha*x(j)*pi*cos(x(j)*pi*((y(j)*(w(i) - 4))/2 - 3)*(x(j) - 1))*(w(i)/2 - 2)*(x(j) - 1);
                -alpha*y(j)*pi*cos(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), -alpha*cos(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))*(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))];
        elseif y(j)>= 0.5 && x(j+1)< 0.5
            J=[ -alpha*cos(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1))*(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1)),   -alpha*x(j)*pi*cos(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1))*(w(i)/2 - 2)*(x(j) - 1);
                alpha*y(j)*pi*cos(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), alpha*cos(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))*(y(j)*pi*((x(j)*(w(i) - 4))/2 - 3) + pi*((x(j)*(w(i) - 4))/2 - 3)*(y(j) - 1))];
        elseif y(j)>= 0.5 && x(j+1)>= 0.5
            J=[ -alpha*cos(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1))*(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1)),   -alpha*x(j)*pi*cos(x(j)*pi*(((w(i) - 4)*(y(j) - 1))/2 + 3)*(x(j) - 1))*(w(i)/2 - 2)*(x(j) - 1);
                -alpha*y(j)*pi*cos(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))*(w(i)/2 - 2)*(y(j) - 1), -alpha*cos(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))*(y(j)*pi*(((w(i) - 4)*(x(j) - 1))/2 + 3) + pi*(((w(i) - 4)*(x(j) - 1))/2 + 3)*(y(j) - 1))];
        end
        
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
plot(w,lamuda(:,1),...
    '-diamond',...
    'color',[144, 44, 252]/255,...
    'MarkerFaceColor',[241, 156, 162]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(w));
hold on ;
plot(w,lamuda(:,2),...
    '-.^',...
    'color',[44, 186, 252]/255,...
    'MarkerFaceColor',[0, 0, 0]/255,...
    'MarkerEdgeColor','none',...
    'MarkerSize',5,...
    'LineWidth',1,...
    'MarkerIndices',1:50:length(w));
y1=line([0,4],[0,0]);
set(y1,'linestyle','--','color','r','LineWidth',1.4);
lgd=legend('\it{LE_1}','\it{LE_2}','Location','SouthEast');
set(lgd,'FontName','Times New Roman','FontSize',17);
set(gca,'LooseInset',get(gca,'TightInset'));
ylabel('LEs','FontName','Times New Roman','FontSize',17);
xlabel('\it{w}','FontName','Times New Roman','FontSize',17);
axis([0 4 -2 1.5])
% saveas(gcf,'NDCM_b.png')

