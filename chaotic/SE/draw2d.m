% ����һЩ���ڻ�ͼ������
[X, Y] = meshgrid(1:1:100,1:1:100);
% Z = -4:6;

% ʹ�� mesh ������������ͼ
mesh(X, Y, SE2D)

% ���� x, y, �� z ������Ŀ̶�
xlim([1,100])
ylim([1,100])
xticks([1,50,100])
yticks([1,50,100])
% zticks(-4:6)
% zlim([0,2.35])
% Ϊ x, y, �� z ��������ӱ�ǩ
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

