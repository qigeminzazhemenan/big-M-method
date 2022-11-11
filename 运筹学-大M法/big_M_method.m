function [] = big_M_method(c,A,symbol_list,max_or_min,b)
% �����η�����׼�����Թ滮����: max cx s.t. Ax=b x>=0
% �������: cΪĿ�꺯��ϵ��, AΪԼ��������ϵ������,symbol_listΪ��������ʽ�Ĳ���ʽ��(1����>,0����=,-1����<)��CbΪԼ�������鳣����
% �������: x���Ž�, z����Ŀ�꺯��ֵ
format rat;
%Ԥ����-������ʼ����---
[m,n] = size(A);  %mԼ����������, n���߱����� 
index_num1 = sum(symbol_list(:)==1)+sum(symbol_list(:)==-1);%index_num1�ɳ�/ʣ���������
index_num2 = sum(symbol_list(:)==1)+sum(symbol_list(:)==0);%index_num2�˹���������

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
label_list = horzcat(index_name_list0,index_name_list1,index_name_list2);%���б�������

label_row = horzcat(["Cb","Xb","b"],label_list,"��");%������ǩ��

c = horzcat(c,repelem(0,index_num1),repelem(-10000*max_or_min,index_num2));%����Ŀ�꺯��ϵ����

A1 = zeros(m,index_num1+index_num2);%A1Լ���������벿��ϵ������
Xb_loc = repelem(0,m);%��ʼ��������λ��(������)
t1 = 1;%t1��������
t2 = index_num1+1;%t2��������
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
A = horzcat(A,A1);%����ϵ������

Xb_loc = Xb_loc+n;%��ʼ��������λ��
Xb = label_list(Xb_loc);%��ʼ������������
Cb = c(Xb_loc);%��ʼ��������Ŀ�꺯���е�ϵ��

L = c-Cb*A;%������

[~,enter_loc] = max(max_or_min*L);
enter_label = label_list(enter_loc);%���������,���������������
xita = (b')./A(:,enter_loc);%������

%������ʼ�����α�
table = horzcat(Cb',Xb',b',A,xita);
table = vertcat(label_row,table);
table = vertcat(horzcat([NaN,NaN,NaN],c,NaN),table);
table = vertcat(table,horzcat([NaN,NaN,NaN],L,NaN))

%ѭ������---    
while max(L*max_or_min)>0

    %����
    xita(xita<=0) = 10000;%��С�ڵ���0ʱ���ÿ���
    [~,out_label] = min(xita);
    out_label = Xb(out_label);%��������������
    out_loc = find(label_list==out_label);%����������λ��
    Xb_loc(Xb_loc==out_loc) = [];
    Xb(Xb==out_label) = [];%����
    Xb_loc(m) = enter_loc;
    Xb(m) = enter_label;%���
    
    %������һ�ֵ�ϵ����������������ֵ
    T = A(:,Xb_loc);%T��������,���ڼ����µ�A��b
    A = T\A;
    b = (T\b')';
    Cb = c(Xb_loc);
    L = c-Cb*A;
    [~,enter_loc] = max(max_or_min*L);
    enter_label = label_list(enter_loc);%���������,���������������
    xita = (b')./A(:,enter_loc);%������

    %������һ�ֵ����α�
    table = horzcat(Cb',Xb',b',A,xita);
    table = vertcat(table,horzcat([NaN,NaN,NaN],L,NaN)) 

end

%������Ž�
if Xb-floor(Xb)~=0.3%�����Ž�
    target_best = Cb*(b')%Ŀ�꺯������ֵ
else%�����Ž�
    target_best = "�����Ž�"%Ŀ�꺯��������ֵ     
end

end

    
    