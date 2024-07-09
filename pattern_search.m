function [x,f_x,count]=pattern_search(x,dim,fobj)
%最后考察参数的具体赋值
count=0;
I=eye(dim); %生成dim阶的单位阵，方便后续单位向量的提取
delta=0.5;%初始步长
alpha=1e-10;%允许误差
a=0.5;%缩减率
b=1.2;%加速系数
j=1;
n=dim;
y=x;%初始点,应该是一个1乘以dim的行向量
%内循环为探测移动，外循环为模式移动
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
        x1=x;%x1相当于算法中的x_k
        x=y;
        y=x+b*(x-x1);
        j=1;%之后进行探测移动
        continue;%转步骤2
    else %移动失败，缩小步长，继续移动
        if delta<=alpha
            break
        else
            delta=delta*a;
            y=x;
            j=1;
%             continue;%转步骤2
        end
    end
end
f_x=fobj(x);




            


