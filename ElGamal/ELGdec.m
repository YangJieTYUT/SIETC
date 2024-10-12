function [x1,re_hash_384] = ELGdec(p,k_E,private_key,y)

re_code_array = repmat('0', length(y), 4);
k_M1 = powermod(k_E,-private_key,p);

x1 = powermod(y .* k_M1,1,p);

for i = 1:length(y)
    re_code_array(i,:) = grayToBinary(dec2bin(x1(i),4));
end

re_hash_384 = transpose(lower(dec2hex(bin2dec(re_code_array))));

end
