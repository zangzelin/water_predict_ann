
clc 
clear all

%% 导入数据
% 计算函数f的雅克比矩阵，是解析式
xlsfile = 'datazzl.txt';                                          %打开文件，输入数据
[data] = xlsread(xlsfile,'a1:e404'); 

%% 数据归一化
 
 mm = mean(data);
 
  for i = 1:5
    traind_s(:,i) = data(:,i) - mm(i);
  end
  
ml(1) = [max(traind_s(:,1))-min(traind_s(:,1))]/2;
ml(2) = [max(traind_s(:,2))-min(traind_s(:,2))]/2;
ml(3) = [max(traind_s(:,3))-min(traind_s(:,3))]/2;
ml(4) = [max(traind_s(:,4))-min(traind_s(:,4))]/2;
ml(5) = [max(traind_s(:,5))-min(traind_s(:,5))]/2;

for i = 1:5
    traind_s(:,i) =  traind_s(:,i)/ml(i);
end


 

% 拟合用数据。参见《数学试验》，p190，例2

data_1=1:404;

obs_1=traind_s(:,1)';

 

% 2. LM算法

% 初始猜测s

a0=10; b0=0.5;

y_init = a0*exp(-b0*data_1);

% 数据个数

Ndata=length(obs_1);

% 参数维数

Nparams=2;

% 迭代最大次数

n_iters=500;

% LM算法的阻尼系数初值

lamda=0.01;

 

% step1: 变量赋值

updateJ=1;

a_est=a0;

b_est=b0;

 

% step2: 迭代

for it=1:n_iters

    if updateJ==1

        % 根据当前估计值，计算雅克比矩阵

        J=zeros(Ndata,Nparams);

        for i=1:length(data_1)

            J(i,:)=[exp(-b_est*data_1(i))                           -a_est*data_1(i)*exp(-b_est*data_1(i))];

        end

        % 根据当前参数，得到函数值

        y_est = a_est*exp(-b_est*data_1);

        % 计算误差

        d=obs_1-y_est;

        % 计算（拟）海塞矩阵

        H=J'*J;

        % 若是第一次迭代，计算误差

        if it==1

            e=dot(d,d);

        end

    end

 

    % 根据阻尼系数lamda混合得到H矩阵

    H_lm=H+(lamda*eye(Nparams,Nparams));

    % 计算步长dp，并根据步长计算新的可能的\参数估计值

    dp=inv(H_lm)*(J'*d(:));

    g = J'*d(:);

    a_lm=a_est+dp(1);

    b_lm=b_est+dp(2);

    % 计算新的可能估计值对应的y和计算残差e

    y_est_lm = a_lm*exp(-b_lm*data_1);
    
    plot(y_est_lm);
    hold on 
    plot(obs_1);
    hold off

    d_lm=obs_1-y_est_lm;

    e_lm=dot(d_lm,d_lm);

    % 根据误差，决定如何更新参数和阻尼系数

    if e_lm<e

        lamda=lamda/10;

        a_est=a_lm;

        b_est=b_lm;

        e=e_lm;

        disp(e);

        updateJ=1;

    else

        updateJ=0;

        lamda=lamda*10;

    end

end

%显示优化的结果

a_est

b_est
