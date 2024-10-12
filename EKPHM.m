function [xn,yn,zn] = EKPHM(key,n)

x(1)= key(1);y(1)=key(2);

p = key(3);
q = key(4);
% r = key(5);
% x(1)= 0.5;y(1)=0.5;
% p = 1;
% q = 1;
% r = 1;
t = 1000;
% n = 100000;
% bait = 0.8;
for j=1:n+t-1
    
    x(j + 1) = mod(exp(pi) * (((p * x(j) * (1 - pi*y(j)) ^ 2)))+ q * y(j) * x(j) ^ 2 + pi * x(j),1);
    y(j + 1) = mod(exp(pi) * (((p * y(j) * (1 - pi*x(j)) ^ 2)))- q * x(j) * y(j) ^ 2 + pi * y(j),1);
%     z(j + 1) = mod(exp(pi) * (((p * z(j) * (1 - pi*x(j)) ^ 2)))- q * x(j) * z(j) ^ 2 + pi * z(j),1);
    
end
xn = x(t+1:end);
yn = y(t+1:end);
% zn = z(t+1:end);
% plot(xn,yn,'.','MarkerSize',1);
end