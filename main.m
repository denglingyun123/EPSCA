clear all
clc
close all
SearchAgents_no=30; % Number of search agents
Max_iteration=1000; % Maximum number of iterations  
Function_name=1; %设定测试函数，1-29.其中大于11的要维度大于10，参考cec2017文档指定的维度
dim=30; %维度设定，维度可供选择范围[2,10,20,30,50,100]，其中Function_name>=11的最低维度设置为10.

lb=-100;%下边界
ub=100;%上边界
fobj = @(x) cec17_func(x',Function_name);

Max_test=30;
for i=1:Max_test
    disp(['第',num2str(i),'次实验']);
    [Best_pos1(i,:),Best_score1(i),SCA_curve(i,:)]=SCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化
    [Best_pos2(i,:),Best_score2(i),EPSCA_curve(i,:)]=EPSCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化 
end
% %结果对比
figure(1)
semilogy(mean(SCA_curve),'color','[1,0.5,0]','linewidth',2.0,'Marker','o','MarkerIndices',1:50:length(mean(SCA_curve)))
hold on
semilogy(mean(EPSCA_curve),'color','[0.19608,0.80392,0.19608]','linewidth',2.0,'Marker','d','MarkerIndices',1:50:length(mean(EPSCA_curve)))


title('Convergence curve of F_{1}')
xlabel('Iteration');
ylabel('Fitness');
axis tight%用 axis tight命令可以让坐标轴调整到紧凑地显示图像或曲线，不留边界的空白
grid off%显示 gca 命令返回的当前坐标区或图的主网格线。主网格线从每个刻度线延伸。
box on %显示坐标区周围的轮廓
legend('SCA','EPSCA')

disp('-------------------------------------------------')
display(['SCA 30次实验最优适应度值(Best) : ', num2str(min(Best_score1))]);
display(['SCA 30次实验最优解对应的平均适应度值(mean) : ', num2str(mean(Best_score1))]);
display(['SCA 30次实验最差适应度值(wrost) : ', num2str(max(Best_score1))]);
display(['SCA 30次实验标准差（std） : ', num2str(std(Best_score1))]);

disp('-------------------------------------------------')
display(['EPSCA 30次实验最优适应度值(Best) : ', num2str(min(Best_score2))]);
display(['EPSCA 30次实验最优解对应的平均适应度值(mean) : ', num2str(mean(Best_score2))]);
display(['EPSCA 30次实验最差适应度值(wrost) : ', num2str(max(Best_score2))]);
display(['EPSCA 30次实验标准差（std） : ', num2str(std(Best_score2))]);
