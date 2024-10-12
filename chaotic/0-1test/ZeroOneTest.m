

function [Kc,Kcreg,KcCorr] =ZeroOneTest(dataSet,cont)

%�������ܣ����0-1����

%�������˵��

%dataSet:����������

%cont:�涨����ֵ,��Χ0-pi��ѡ

%------------------------

%�������˵��

%Kc������Mc����

%Kcreg����������D�õ�

%------------------------

%��������

%dataset = logisticmap(0.3,3.97,4999);

%cont = 2;

% [Kc,Kcreg,KcCorr] =ZeroOneTest(dataSet,cont);

%�ο����ף�The 0-1 Test for Chaos: A Review

%���㿪ʼ

Ndata = size(dataSet,1);

p(1)=dataSet(1)*cos(cont);

s(1)=dataSet(1)*sin(cont);

for n = 1:Ndata
    
    
    p(n+1)=p(n)+dataSet(n)*cos(n*cont);
    
    
    s(n+1)=s(n)+dataSet(n)*sin(n*cont);
    
end

% %����p��s�Ĺ켣ͼ

% figure

% plot(p,s);

% xlabel('p_c');ylabel('q_c');

% %��ͼ������

%�������λ��

Numn = Ndata/10;%���鲻�������ݼ���С��ʮ��ʮһ

meanpower = mean(dataSet,1)^2;

for n = 1:Numn
    
    sumMean = 0;
    
    for j = 1:Ndata-Numn%�˴���ʽ��N�Ĵ�Сȡʣ���ֵ
        
        
        sumMean = sumMean + (p(j+n)-p(j))^2 + (s(j+n)-s(j))^2;
        
    end
    
    Mc(n) = 1/Ndata*sumMean;
    
    Vosc(n) = meanpower*((1-cos(n*cont))/(1-cos(cont)));
    
end

D = Mc -Vosc;

x = 1:Numn;

KcCorr = corr(x',D');

% %��ͼ������Mc��D�ı仯

% figure

% n = 1:Numn;

% plot(n,Mc);

% hold on

% plot(n,D);

% xlabel('n');ylabel('Mc,D');

% %������ͼ

Dreg = D(Numn) + 1.1*min(abs(D));

Kc = log(Mc(Numn))/log(Numn);

Kcreg = log(Dreg)/log(Numn);

end
