function gray = binaryToGray(binary)
    % ���������ַ���ת��Ϊ��������
    binaryNum = binary - '0';
    
    % ��ʼ����������������
    grayNum = zeros(1, length(binary));

    % ��һλ����
    grayNum(1) = binaryNum(1);

    % ��ʣ���λ����������
    for i = 2:length(binaryNum)
        grayNum(i) = xor(binaryNum(i - 1), binaryNum(i));
    end

    % ����������ת�����ַ���
    gray = char(grayNum + '0');
end