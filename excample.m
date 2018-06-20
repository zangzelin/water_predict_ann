
clc 
clear all

%% ��������
% ���㺯��f���ſ˱Ⱦ����ǽ���ʽ
xlsfile = 'datazzl.txt';                                          %���ļ�����������
[data] = xlsread(xlsfile,'a1:e404'); 

%% ���ݹ�һ��
 
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


 

% ��������ݡ��μ�����ѧ���顷��p190����2

data_1=1:404;

obs_1=traind_s(:,1)';

 

% 2. LM�㷨

% ��ʼ�²�s

a0=10; b0=0.5;

y_init = a0*exp(-b0*data_1);

% ���ݸ���

Ndata=length(obs_1);

% ����ά��

Nparams=2;

% ����������

n_iters=500;

% LM�㷨������ϵ����ֵ

lamda=0.01;

 

% step1: ������ֵ

updateJ=1;

a_est=a0;

b_est=b0;

 

% step2: ����

for it=1:n_iters

    if updateJ==1

        % ���ݵ�ǰ����ֵ�������ſ˱Ⱦ���

        J=zeros(Ndata,Nparams);

        for i=1:length(data_1)

            J(i,:)=[exp(-b_est*data_1(i))                           -a_est*data_1(i)*exp(-b_est*data_1(i))];

        end

        % ���ݵ�ǰ�������õ�����ֵ

        y_est = a_est*exp(-b_est*data_1);

        % �������

        d=obs_1-y_est;

        % ���㣨�⣩��������

        H=J'*J;

        % ���ǵ�һ�ε������������

        if it==1

            e=dot(d,d);

        end

    end

 

    % ��������ϵ��lamda��ϵõ�H����

    H_lm=H+(lamda*eye(Nparams,Nparams));

    % ���㲽��dp�������ݲ��������µĿ��ܵ�\��������ֵ

    dp=inv(H_lm)*(J'*d(:));

    g = J'*d(:);

    a_lm=a_est+dp(1);

    b_lm=b_est+dp(2);

    % �����µĿ��ܹ���ֵ��Ӧ��y�ͼ���в�e

    y_est_lm = a_lm*exp(-b_lm*data_1);
    
    plot(y_est_lm);
    hold on 
    plot(obs_1);
    hold off

    d_lm=obs_1-y_est_lm;

    e_lm=dot(d_lm,d_lm);

    % ������������θ��²���������ϵ��

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

%��ʾ�Ż��Ľ��

a_est

b_est
