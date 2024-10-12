function h = hash(P,meth)

P=P(:);
% convert strings and logicals into uint8 format
if ischar(P) || islogical(P)
    P=uint8(P);
else % convert everything else into uint8 format without loss of data
    P=typecast(P,'uint8');
end

% verify hash method, with some syntactical forgiveness:
meth=upper(meth);
switch meth
    case 'SHA1'
        meth='SHA-1';
    case 'SHA256'
        meth='SHA-256';
    case 'SHA384'
        meth='SHA-384';
    case 'SHA512'
        meth='SHA-512';
    otherwise
end
algs={'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if isempty(strcmp(algs,meth))
    error(['Hash algorithm must be ' ...
        'MD2, MD5, SHA-1, SHA-256, SHA-384, or SHA-512']);
end

% create hash
x=java.security.MessageDigest.getInstance(meth);
x.update(P);
h=typecast(x.digest,'uint8');
h=dec2hex(h)';
if(size(h,1))==1 % remote possibility: all hash bytes  128, so pad:
    h=[repmat('0',[1 size(h,2)]);h];
end
h=lower(h(:)');
clear x
return


