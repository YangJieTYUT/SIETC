% 创建一些用于绘图的数据
[X, Y] = meshgrid(1:1:100,1:1:100);
% Z = -4:6;

% 使用 mesh 函数创建网格图
mesh(X, Y, SE2D)

% 调整 x, y, 和 z 坐标轴的刻度
xlim([1,100])
ylim([1,100])
xticks([1,50,100])
yticks([1,50,100])
% zticks(-4:6)
% zlim([0,2.35])
% 为 x, y, 和 z 坐标轴添加标签
% xlabel('\it{\mu}')
% ylabel('\it{k}')
xlabel('\it{p}')
ylabel('\it{q}')
% xlabel('\it{c}')
% ylabel('\it{a}')
% xlabel('\it{a}')
% ylabel('\it{\mu}')
% 3
% xlabel('\it{a}')
% ylabel('\it{k}')
zlabel('SE')

