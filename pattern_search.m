function [x,f_x,count]=pattern_search(x,dim,fobj)
%��󿼲�����ľ��帳ֵ
count=0;
I=eye(dim); %����dim�׵ĵ�λ�󣬷��������λ��������ȡ
delta=0.5;%��ʼ����
alpha=1e-10;%�������
a=0.5;%������
b=1.2;%����ϵ��
j=1;
n=dim;
y=x;%��ʼ��,Ӧ����һ��1����dim��������
%��ѭ��Ϊ̽���ƶ�����ѭ��Ϊģʽ�ƶ�
while delta>alpha
    while j<=n
        f1=fobj(y+delta*I(j,:));
        count=count+1;
        if f1<fobj(y)
            y=y+delta*I(j,:);
            j=j+1;
            continue;
        end
        
        f2=fobj(y-delta*I(j,:));
        count=count+1;
        if f2<fobj(y)
            y=y-delta*I(j,:);
            j=j+1;
            continue;
        else
            j=j+1;
        end
    end
    
    f3=fobj(y);
    count=count+1;
    if f3<fobj(x)
        x1=x;%x1�൱���㷨�е�x_k
        x=y;
        y=x+b*(x-x1);
        j=1;%֮�����̽���ƶ�
        continue;%ת����2
    else %�ƶ�ʧ�ܣ���С�����������ƶ�
        if delta<=alpha
            break
        else
            delta=delta*a;
            y=x;
            j=1;
%             continue;%ת����2
        end
    end
end
f_x=fobj(x);




            


