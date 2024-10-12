function binary = grayToBinary(gray)
    % ���������ַ���ת��Ϊ��������
    grayNum = gray - '0';
    
    % ��ʼ������������������
    binaryNum = zeros(1, length(gray));

    % ��һλ����
    binaryNum(1) = grayNum(1);

    % ��ʣ���λ������λ������
    for i = 2:length(grayNum)
        binaryNum(i) = xor(binaryNum(i - 1), grayNum(i));
    end

    % ����������ת�����ַ���
    binary = char(binaryNum + '0');
end