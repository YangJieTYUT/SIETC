% clear;clc;close all;

n=1000;
x1(1)=0.6;y1(1)=0.7;
x2(1)=0.6;y2(1)=0.7;
x3(1)=0.6;y3(1)=0.7;
x4(1)=0.6;y4(1)=0.7;

m=5;t=1;
SE1=[];C1=[];SE2=[];C2=[];SE3=[];C3=[];SE4=[];C4=[];

% new EPHM system

for p=0:1:100
    for q = 0:1:100
    for i=1:n-1 
        	x4(i + 1) = mod(exp(pi) * (((p * x4(i) * (1 - pi*y4(i)) ^ 2)))+ q * y4(i) * x4(i) ^ 2 + pi * x4(i),1);
            y4(i + 1) = mod(exp(pi) * (((p * y4(i) * (1 - pi*x4(i)) ^ 2)))- q * x4(i) * y4(i) ^ 2 + pi * y4(i),1);
    end
   SE4 = [SE4,SampEn(x4,2,0.2*std(x4))];
%    C4=[C4;p];
    end
end

% 2D-SLTCM
for w=0.001:0.01:4
    alpha = 1;
    for i=1:n-1
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
    end
    SE3 = [SE3,SampEn(x3,2,0.2*std(x3))];
    C3=[C3;w];
end

% new discrete chaotic map
for b=0.001:0.01:4
    a = 0.8;k = 0.5;
    for i=1:n-1
        x2(i+1)= (1-a * x2(i)^2 + y2(i)) * sin(k/x2(i));
        y2(i+1)= b* x2(i)* sin(k/x2(i));
    end
    SE2 = [SE2,SampEn(x2,2,0.2*std(x2))];
    C2=[C2;b];
end

% improved 2D-LASM
for mu=0.001:0.01:4
    rho = 4;
    for i=1:n-1 
        x1(i+1)=sin(pi * mu * (y1(i) + rho - x1(i)) * x1(i) * (1 - x1(i)));
        y1(i+1)=sin(pi * mu * (x1(i+1) + rho - y1(i)) * y1(i) * (1 - y1(i)));
    end
   SE1 = [SE1,SampEn(x1,2,0.2*std(x1))];
   C1=[C1;mu];
end


%%
close all

p1=plot(C1,SE1,'--b*','linewidth',1,'MarkerSize',5,'MarkerIndices',1:15:length(SE1));hold on
p2=plot(C2,SE2,'k-.pentagram','linewidth',2,'MarkerFaceColor','m','MarkerSize',5,'MarkerIndices',1:15:length(SE1));hold on
p3=plot(C3,SE3,':cdiamond','linewidth',2,'MarkerFaceColor','m','MarkerSize',5,'MarkerIndices',1:15:length(SE1));hold on
p4=plot(C4,SE4,'-rhexagram','linewidth',2,'MarkerFaceColor','r','MarkerSize',5,'MarkerIndices',1:15:length(SE1));

set(p2,'color',[142 207 201]/255,'LineWidth',1);
set(p3,'color',[255 190 122]/255,'LineWidth',1);
set(p4,'color',[130 176 210]/255,'LineWidth',1);

% set(p5,'color',[250 127 111]/255,'LineWidth',1.5);



set(gca,'FontName','Times New Roman');
set(gca,'LooseInset',get(gca,'TightInset'),'linewidth',1);
lgd=legend('2D-ILASM','2D-NDCM','2D-SLTCM','2D-EPHM','location','southeast'); 


set(lgd,'FontSize',10);
set(gca,'LooseInset',get(gca,'TightInset'));
% set(gca,'YLim',[0 1],'FontSize',17);set(gca,'YTick',0:0.2:1);set(gca,'YTickLabel',0:0.2:1);
xlabel('\mu,\it{b},\it{w},\it{p}');ylabel('SEs');

% legend boxoff;

% saveas(gcf,'PEs.png')
% grid on;











