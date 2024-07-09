clear all
clc
close all
SearchAgents_no=30; % Number of search agents
Max_iteration=1000; % Maximum number of iterations  
Function_name=1; %�趨���Ժ�����1-29.���д���11��Ҫά�ȴ���10���ο�cec2017�ĵ�ָ����ά��
dim=30; %ά���趨��ά�ȿɹ�ѡ��Χ[2,10,20,30,50,100]������Function_name>=11�����ά������Ϊ10.

lb=-100;%�±߽�
ub=100;%�ϱ߽�
fobj = @(x) cec17_func(x',Function_name);

Max_test=30;
for i=1:Max_test
    disp(['��',num2str(i),'��ʵ��']);
    [Best_pos1(i,:),Best_score1(i),SCA_curve(i,:)]=SCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %��ʼ�Ż�
    [Best_pos2(i,:),Best_score2(i),EPSCA_curve(i,:)]=EPSCA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %��ʼ�Ż� 
end
% %����Ա�
figure(1)
semilogy(mean(SCA_curve),'color','[1,0.5,0]','linewidth',2.0,'Marker','o','MarkerIndices',1:50:length(mean(SCA_curve)))
hold on
semilogy(mean(EPSCA_curve),'color','[0.19608,0.80392,0.19608]','linewidth',2.0,'Marker','d','MarkerIndices',1:50:length(mean(EPSCA_curve)))


title('Convergence curve of F_{1}')
xlabel('Iteration');
ylabel('Fitness');
axis tight%�� axis tight�����������������������յ���ʾͼ������ߣ������߽�Ŀհ�
grid off%��ʾ gca ����صĵ�ǰ��������ͼ���������ߡ��������ߴ�ÿ���̶������졣
box on %��ʾ��������Χ������
legend('SCA','EPSCA')

disp('-------------------------------------------------')
display(['SCA 30��ʵ��������Ӧ��ֵ(Best) : ', num2str(min(Best_score1))]);
display(['SCA 30��ʵ�����Ž��Ӧ��ƽ����Ӧ��ֵ(mean) : ', num2str(mean(Best_score1))]);
display(['SCA 30��ʵ�������Ӧ��ֵ(wrost) : ', num2str(max(Best_score1))]);
display(['SCA 30��ʵ���׼�std�� : ', num2str(std(Best_score1))]);

disp('-------------------------------------------------')
display(['EPSCA 30��ʵ��������Ӧ��ֵ(Best) : ', num2str(min(Best_score2))]);
display(['EPSCA 30��ʵ�����Ž��Ӧ��ƽ����Ӧ��ֵ(mean) : ', num2str(mean(Best_score2))]);
display(['EPSCA 30��ʵ�������Ӧ��ֵ(wrost) : ', num2str(max(Best_score2))]);
display(['EPSCA 30��ʵ���׼�std�� : ', num2str(std(Best_score2))]);
