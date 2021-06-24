%% I. 清空环境变量
clear all
clc

%% II. 导入数据
citys = [716180.123,3530667.308,70.011;
716180.123,3530668.20000214,60.412;
716180.123,3530668.308,50.819;
716180.123,3530669.20000214,43.360;
716180.123,3530670.308,40.849;
716180.123,3530671.20000214,64.954;
716180.123,3530672.308,73.021;
716177.123,3530694.20000214,64.747;
716180.123,3530718.225,73.021;
716180.123,3530719.308,65.015;
716180.123,3530720.654,55.729];
%citys = [2 5 91;3 6 85;4 7 55;5 8 45;7 9 25;9 10 70;...
 %   11 11 92; 20 12 95;29 13 94;31 14 71;35 15 50;0 0 0];

V_LAND =10; %地速，起飞前为无人机设定的速度
V_WIND =5; %风速大小
POWER = 166;
beta_wind = 5*pi/4;%风速的夹角

%% III. 计算城市间相互距离
n = size(citys,1);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            time_direct = sqrt(sum((citys(i,:) - citys(j,:)).^2))/V_LAND;
        		D(i,j) = calculate_V(citys(j,3),citys(j,2),citys(j,1),citys(i,3),citys(i,2), citys(i,1))*POWER*time_direct ;
                D(j,i) = calculate_V(citys(i,3),citys(i,2), citys(i,1), citys(j,3),citys(j,2),citys(j,1))*POWER*time_direct ;
        else
            D(i,j) = 1e-4;      
        end
    end    
end

%% IV. 初始化参数
m = 100;                             % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
rho = 0.2;                           % 信息素挥发因子
Q = 10;                              % 常系数
Eta = 1./D;                          % 启发函数
Tau = ones(n,n);                     % 信息素矩阵
Table = zeros(m,n);                  % 路径记录表
iter = 1;                            % 迭代次数初值
iter_max = 500;                      % 最大迭代次数 
Route_best = zeros(iter_max,n);      % 各代最佳路径       
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度
Limit_iter = 0;                      % 程序收敛时迭代次数

%% V. 迭代寻找最佳路径
while iter <= iter_max
     % 随机产生各个蚂蚁的起点城市
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(n);
          start(i) = temp(1);%50个蚂蚁的出发城市
          %start(i) = 12;
      end
      Table(:,1) = start; 
      citys_index = 1:n;
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个城市路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % 已访问的城市集合(禁忌表)
             allow_index = ~ismember(citys_index,tabu);%ismember检查第一个矩阵中元素是否为第二个矩阵中的值
             allow = citys_index(allow_index);  % 待访问的城市集合
             P = allow;
             % 计算城市间转移概率
             for k = 1:length(allow)
                 P(k) = (Tau(tabu(end),allow(k))^alpha) * (Eta(tabu(end),allow(k))^beta);
             end
             P = P/sum(P);
             % 轮盘赌法选择下一个访问城市
             Pc = cumsum(P);     
            target_index = find(Pc >= rand); 
            target = allow(target_index(1));
            Table(i,j) = target;
         end
      end
      % 计算各个蚂蚁的路径距离
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % 计算最短路径距离及平均距离
      if iter == 1
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);
          Limit_iter = 1;
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
              Limit_iter = iter;
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个城市计算
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
    % 迭代次数加1，清空路径记录表
    %Rlength(iter) = min_Length;
    iter = iter + 1;
    Table = zeros(m,n);
end

%% VI. 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);

%% VII. 绘图
figure(1)
plot3([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],...
     [citys(Shortest_Route,3);citys(Shortest_Route(1),3),],'o-');
 
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),citys(i,3),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),citys(Shortest_Route(1),3),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),citys(Shortest_Route(end),3),'       终点');
set(gca, 'xtick', 0 : 5 : 40);
set(gca, 'ytick',0 : 2 : 16, 'ydir','reverse' );
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
zlabel('城市位置高程坐标')
title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('各代最短距离与平均距离对比')
figure(3)
plot(1:iter_max,Length_best,'b')
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')