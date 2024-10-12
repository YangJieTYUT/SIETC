function binary = grayToBinary(gray)
    % 将格雷码字符串转换为数字数组
    grayNum = gray - '0';
    
    % 初始化二进制码数字数组
    binaryNum = zeros(1, length(gray));

    % 第一位不变
    binaryNum(1) = grayNum(1);

    % 对剩余的位进行逐位异或操作
    for i = 2:length(grayNum)
        binaryNum(i) = xor(binaryNum(i - 1), grayNum(i));
    end

    % 将数字数组转换回字符串
    binary = char(binaryNum + '0');
end