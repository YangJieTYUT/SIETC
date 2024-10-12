function [key1,key2,key3] = key_gen(sequence,char)
    hash_value = hash(sequence,char);
    hash_512 = hex2dec(hash_value');
    x0 = mod((bin2dec(reshape(dec2bin(bitxor(hash_512(1:10),hash_512(11:20)),4)',1,[]))+bin2dec(reshape(dec2bin(bitxor(hash_512(21:25),hash_512(26:30)))',1,[])))/1e11,1);
    y0 = mod((bin2dec(reshape(dec2bin(bitxor(hash_512(31:40),hash_512(41:50)),4)',1,[]))+bin2dec(reshape(dec2bin(bitxor(hash_512(51:55),hash_512(56:60)))',1,[])))/1e11,1);
    p0 = (bin2dec(reshape(dec2bin(hash_512(61:70))',1,[]))+0.01* bin2dec(reshape(dec2bin(hash_512(71:80))',1,[])))/1e11;
    q0 = (bin2dec(reshape(dec2bin(hash_512(81:90))',1,[]))+0.01* bin2dec(reshape(dec2bin(hash_512(91:100))',1,[])))/1e11;
    alpha1 = sum(hash_512(101:114))/max(hash_512(115:128));
    alpha2 = sum(hash_512(115:128))/max(hash_512(101:114));
    
    key1 = [x0,y0,p0,q0];
    key2 = [x0,y0,p0+alpha1,q0 + alpha2];
    key3 = [mod(x0+alpha1,1),mod(y0 + alpha2,1),p0,q0];
end
