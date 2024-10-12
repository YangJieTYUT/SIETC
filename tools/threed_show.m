function threed_show( image_data )
% ��3D����ת��Ϊ��������Ƭ
slices = num2cell(image_data, [1 2]);

% ʹ��arrayfun����ÿ����Ƭ
rotated_slices = arrayfun(@(x) imrotate(x{1}, 180, 'bilinear', 'crop'), slices, 'UniformOutput', false);

% ����ת�����Ƭ��cell����ת��3D����
image_data = cat(3, rotated_slices{:});
%
figure;

% ���� img ��һ�� MxNxP ����άͼ�����
[M, N, P] = size(image_data);

% ��������������Ӧ3D����
axis([1 M 1 N 1 P]);
% axis off
hold on;

for k = 1:P
    % ����ÿ����Ƭ������ʹ��һ�� surf ��������ʾ
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

