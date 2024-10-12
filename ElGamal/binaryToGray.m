function gray = binaryToGray(binary)
    % 将二进制字符串转换为数字数组
    binaryNum = binary - '0';
    
    % 初始化格雷码数字数组
    grayNum = zeros(1, length(binary));

    % 第一位不变
    grayNum(1) = binaryNum(1);

    % 对剩余的位进行异或操作
    for i = 2:length(binaryNum)
        grayNum(i) = xor(binaryNum(i - 1), binaryNum(i));
    end

    % 将数字数组转换回字符串
    gray = char(grayNum + '0');
end