function [public_key,k_E,y] = ELGenc(p,private_key,hash_384)
%������ϢΪ˽Կ����Կ��ԭʼ���������Լ�Ҫ���ܵ���Ϣhashֵ

hash_10 = zeros(1,length(hash_384));
inix = zeros(1,length(hash_384));
grayCodeArray = repmat('0', length(hash_384), 4);


for i = 1:length(hash_384)
    hash_10(i) = hex2dec(hash_384(i));
    grayCodeArray(i,:) = binaryToGray(dec2bin(hash_10(i),4));
    inix(i) = bin2dec(grayCodeArray(i,:));
end

alpha = primitiveroot(p);

public_key = powermod(alpha^private_key,1,p);

i  = 2;
k_E = powermod(alpha^i,1,p);

y = powermod(inix .* public_key^i,1,p);
end