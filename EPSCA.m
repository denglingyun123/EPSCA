%  A Sine Cosine Algorithm guided by elite pool (EPSCA)                                                                                                                                                                 
%  Developed in MATLAB R2018b                                                                                                                                                                    
%  Authors: Lingyun Deng and Sanyang Liu                                                      
                                                                                                      
%  Main paper:                                                                                        
%   SCA: A Sine Cosine Algorithm guided by elite pool strategy for global optimization
%  Applied soft computing, DOI: https://doi.org/10.1016/j.asoc.2024.111946

function [Destination_position,Destination_fitness,Convergence_curve]=EPSCA(N,Max_iter,lb,ub,dim,fobj)
% tic
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
%Initialize the set of random solutions
X=initialization_SSA(N,dim,ub,lb);
% Div=[];%用于存储每一代的多样性值
% epl=[];%用于存储每一代的探索率
% ept=[];%用于存储每一代的开发率
% Div_temp=[];
% for j=1:dim
%     sum1=0;
%     for i=1:N
%         sum1=sum1+abs(median(X(:,j))-X(i,j));
%     end
%     Div_temp(j)=sum1/N;
% end
% Div(1)=sum(Div_temp)/dim;
% epl(1)=(Div(1)/max(Div))*100;
% ept(1)=(abs(Div(1)-max(Div))/max(Div))*100;

Destination_position=zeros(1,dim);
Destination_fitness=inf;
Objective_values = zeros(1,size(X,1));

Convergence_curve=[];
N1=N/2;
Elite_pool=[];
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end

[~,idx1]=sort(Objective_values);
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
sum1=0;
for i=1:N1
    sum1=sum1+X(idx1(i),:);
end
half_best_mean=sum1/N1;
Elite_pool(1,:)=Destination_position;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;


Convergence_curve(1) = Destination_fitness;
%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while l<=Max_iter
    a = 2;
    r1=a-l*((a)/Max_iter); % r1 decreases linearly from a to 0
    RB=randn(N,dim);          %Brownian random number vector
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        k1=randperm(4,1);
        for j=1:size(X,2) % in j-th dimension
            r2=(2*pi)*rand();
            r3=rand();
          %% Elite pool strategy and Brownian motion
            if r3<0.5
                X(i,j)= Elite_pool(k1,j)+RB(i,j)*(r1*sin(r2)*abs(Destination_position(j)-X(i,j)));
            else
                X(i,j)= Elite_pool(k1,j)+RB(i,j)*(r1*cos(r2)*abs(Destination_position(j)-X(i,j)));
            end
        end
    end
    
    % Check if solutions go outside the search spaceand bring them back
    for i=1:size(X,1)
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    %计算探索率与开发率
%    for j=1:dim
%        sum1=0;
%        for i=1:N
%            sum1=sum1+abs(median(X(:,j))-X(i,j));
%        end
%        Div_temp(j)=sum1/N;
%    end
%    Div(l)=sum(Div_temp)/dim;
%    epl(l)=(Div(l)/max(Div))*100;
%    ept(l)=(abs(Div(l)-max(Div))/max(Div))*100;
    
    
    %% pattern search
    if rem(l,90)==0
        [GBestX_temp,~]=pattern_search(Destination_position,dim,fobj);
        %边界控制
       for j = 1: dim
           if(GBestX_temp(j)>ub(j))
               GBestX_temp(j) =ub(j);
           end
           if(GBestX_temp(j)<lb(j))
               GBestX_temp(j) =lb(j);
           end
       end
       %评价
       f=fobj(GBestX_temp);
       if f<Destination_fitness
           Destination_position=GBestX_temp;
           Destination_fitness=f;
       end
    end

   %% An efficient mutation operator
   X_new1=lb+rand*(ub-lb);
   f_new1=fobj(X_new1);
   k=rand;
   X_new2=k*Destination_position+(1-k)*X_new1;
   f_new2=fobj(X_new2);
   % Improve the quality of the current best solution by the mutation operator
   if f_new1<Destination_fitness && f_new2<Destination_fitness && f_new1<=f_new2
       Destination_position=X_new1;
       Destination_fitness=f_new1;
   end
   if f_new1<Destination_fitness && f_new2<Destination_fitness && f_new2<=f_new1
       Destination_position=X_new2;
       Destination_fitness=f_new2;
   end
   if f_new1<Destination_fitness && f_new2>=Destination_fitness
       Destination_position=X_new1;
       Destination_fitness=f_new1;
   end
   if f_new1>=Destination_fitness && f_new2<Destination_fitness
       Destination_position=X_new2;
       Destination_fitness=f_new2;
   end
   
    %% Update the elite pool
    [~,idx1]=sort(Objective_values);
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    sum1=0;
    for i=1:N1
        sum1=sum1+X(idx1(i),:);
    end
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Destination_position;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;

    Convergence_curve(l)=Destination_fitness;
    l=l+1;
end
% toc
% disp(['The runtime of EPSCA: ',num2str(toc)]);
