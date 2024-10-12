function threed_show( image_data )
% 将3D矩阵转换为单独的切片
slices = num2cell(image_data, [1 2]);

% 使用arrayfun处理每个切片
rotated_slices = arrayfun(@(x) imrotate(x{1}, 180, 'bilinear', 'crop'), slices, 'UniformOutput', false);

% 将旋转后的切片从cell数组转回3D矩阵
image_data = cat(3, rotated_slices{:});
%
figure;

% 假设 img 是一个 MxNxP 的三维图像矩阵
[M, N, P] = size(image_data);

% 调整坐标轴以适应3D数据
axis([1 M 1 N 1 P]);
% axis off
hold on;

for k = 1:P
    % 对于每个切片，我们使用一个 surf 对象来显示
    x = [1 N; 1 N];
    y = [1 1; M M];
    z = k * ones(2);
    
    surf(z,x, y, 'CData', image_data(:,:,P-k+1), 'FaceColor', 'texturemap', 'EdgeColor', 'none');
end

view(45, 30);
colormap gray;
axis tight;
if M== 256
    zlim([1,256])
    zticks([1,50,100,150,200,256])
    ylim([1,256])
    yticks([1,50,100,150,200,256])
elseif M== 512
    zlim([1,512])
    zticks([1,100,200,300,400,512])
    ylim([1,512])
    yticks([1,100,200,300,400,512])
end
end

