function pe = pec(y,m,t)

% m=5;t=1;
ly = length(y);
permlist = perms(1:m);%ȫ����m��
c(1:length(permlist))=0;%
    
 for j=1:ly-t*(m-1)%K=n-(m-1)*t ��k���ع�����
     [a,iv]=sort(y(j:t:j+t*(m-1)));%������ֵ��С�������� iv����ֵ aΪ��������ֵ
     for jj=1:length(permlist)%m��
         if (abs(permlist(jj,:)-iv))==0% iv�ɿ����õ��ķ�������
             c(jj) = c(jj) + 1 ;%�ۼ�ÿ�����ж�Ӧm���г��ֵĴ���
         end
     end
 end
 
c=c(find(c~=0));%�ҳ�����Ԫ��������ֵ��Ӧ����ֵ
p = c/sum(c);%����ÿһ�ַ������г��ֵĸ���
pe = -sum(p .* log(p));%shannon�ص���ʽ
pe=pe/log(factorial(m));
% P_E=[P_E;pe/log(factorial(m));]
% C=[C;st];
end