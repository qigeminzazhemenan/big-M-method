function [] = big_M_method(c,A,symbol_list,max_or_min,b)
% 单纯形法求解标准形线性规划问题: max cx s.t. Ax=b x>=0
% 输入参数: c为目标函数系数, A为约束方程组系数矩阵,symbol_list为各个不等式的不等式号(1代表>,0代表=,-1代表<)，Cb为约束方程组常数项
% 输出参数: x最优解, z最优目标函数值
format rat;
%预处理-创建初始矩阵---
[m,n] = size(A);  %m约束条件个数, n决策变量数 
index_num1 = sum(symbol_list(:)==1)+sum(symbol_list(:)==-1);%index_num1松弛/剩余变量个数
index_num2 = sum(symbol_list(:)==1)+sum(symbol_list(:)==0);%index_num2人工变量个数

index_name_list0 = repelem(0.1,n);
index_name_list1 = repelem(0.2,index_num1);
index_name_list2 = repelem(0.3,index_num2);
for i=1:n
    index_name_list0(i) = index_name_list0(i)+i;
end
for i=1:index_num1
    index_name_list1(i) = index_name_list1(i)+i;
end
for i=1:index_num2
    index_name_list2(i) = index_name_list2(i)+i;
end
label_list = horzcat(index_name_list0,index_name_list1,index_name_list2);%所有变量名称

label_row = horzcat(["Cb","Xb","b"],label_list,"θ");%创建标签行

c = horzcat(c,repelem(0,index_num1),repelem(-10000*max_or_min,index_num2));%扩充目标函数系数行

A1 = zeros(m,index_num1+index_num2);%A1约束方程组后半部分系数矩阵
Xb_loc = repelem(0,m);%初始基变量的位置(非最终)
t1 = 1;%t1辅助索引
t2 = index_num1+1;%t2辅助索引
for i =1:m
    if symbol_list(i)==1
        A1(i,t1) = -1;
        A1(i,t2) = 1;
        Xb_loc(i) = t2;
        t1 = t1+1;
        t2 = t2+1;
    elseif symbol_list(i)==0
        A1(i,t2) = 1;
        Xb_loc(i) = t2;
        t2 = t2+1;
    elseif symbol_list(i)==-1
        A1(i,t1) = 1;
        Xb_loc(i) = t1;
        t1 = t1+1;
    end
end
A = horzcat(A,A1);%扩充系数矩阵

Xb_loc = Xb_loc+n;%初始基变量的位置
Xb = label_list(Xb_loc);%初始基变量的名称
Cb = c(Xb_loc);%初始基变量在目标函数中的系数

L = c-Cb*A;%检验数

[~,enter_loc] = max(max_or_min*L);
enter_label = label_list(enter_loc);%检验数最大,即入基变量的名称
xita = (b')./A(:,enter_loc);%西塔列

%构建初始单纯形表
table = horzcat(Cb',Xb',b',A,xita);
table = vertcat(label_row,table);
table = vertcat(horzcat([NaN,NaN,NaN],c,NaN),table);
table = vertcat(table,horzcat([NaN,NaN,NaN],L,NaN))

%循环换基---    
while max(L*max_or_min)>0

    %换基
    xita(xita<=0) = 10000;%θ小于等于0时不用考虑
    [~,out_label] = min(xita);
    out_label = Xb(out_label);%出基变量的名称
    out_loc = find(label_list==out_label);%出基变量的位置
    Xb_loc(Xb_loc==out_loc) = [];
    Xb(Xb==out_label) = [];%出基
    Xb_loc(m) = enter_loc;
    Xb(m) = enter_label;%入基
    
    %计算新一轮的系数，检验数，西塔值
    T = A(:,Xb_loc);%T辅助矩阵,用于计算新的A与b
    A = T\A;
    b = (T\b')';
    Cb = c(Xb_loc);
    L = c-Cb*A;
    [~,enter_loc] = max(max_or_min*L);
    enter_label = label_list(enter_loc);%检验数最大,即入基变量的名称
    xita = (b')./A(:,enter_loc);%西塔列

    %构建新一轮单纯形表
    table = horzcat(Cb',Xb',b',A,xita);
    table = vertcat(table,horzcat([NaN,NaN,NaN],L,NaN)) 

end

%输出最优解
if Xb-floor(Xb)~=0.3%有最优解
    target_best = Cb*(b')%目标函数最优值
else%无最优解
    target_best = "无最优解"%目标函数无最优值     
end

end

    
    