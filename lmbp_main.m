% bp�㷨
%���㷨ʹ��BP�������ˮ���������Ԥ��
%����Ϊ�Ǹ��������ڵ�ˮ������
%���Ϊ��ʮһ�����ڵ�ˮ������
%ѵ��ʹ�õ����ݼ���404��ˮ������
% ���ߣ������ʿ
% ���ڣ�2016/5/1
% �汾��v1

%% ״̬��� 
clc 
clear all

%% ��������

xlsfile = 'datazzl.txt';                                          %���ļ�����������
[data] = xlsread(xlsfile,'a1:e404'); 

%% ���ݹ�һ��

  for i = 1:4
    
    ml(i) = [max(data(:,i))-min(data(:,i))]/2;
    data(:,i) =  data(:,i)/ml(i);
     
    mm(i) = mean(data(:,i));
    traind_s(:,i) = data(:,i) - mm(i);
    
  end
   % traind_s(:,5) = traind_s(:,5) +0.8;

%%  �����������

eb = 0.0001 ;                                                  		%�������
eta = 0.15;                                                    		%ѧϰ��
mc =0.00000;                                               			%����
mcb = 0;

maxiter = 404 ;                                                     %������

iteration = 0;                                               		%��һ�ε���
ntrainnum =404;
inbox = 10;
%% ���繹��

net.nin = 40;                                                         
net.nhidden = 39;
net.nout = 4;                                                      	%���������
dWEXold = 0;
dwexold = 0;

w = 2*(rand(net.nin,net.nhidden)-1/2);
b = 2*(rand(net.nin,1)-1/2);
net.w1 = [w,b];                                                   	%�����Ȩֵ
uk = 1;

W= 2*(rand(net.nhidden+2,net.nout)-1/2);

net.w2 =[W]     ;                                          			%������Ȩֵ
%% ��������

errrec = zeros(1,4);                                           		%������
outrec = zeros(4,maxiter);                                			%����¼
errsum(1,:) = zeros(1,4);
NET = [ ];
index = 1;

%% ����½�
for bianshu =1:150
    
    for i = 1 : maxiter-10
        %% ���ļ���
        
        %�������н�ȡ���ε�����������
        sampinex1 = traind_s(i : i + inbox - 1,:);
        for x = 1:10
        sampinex2(1,4*(x-1)+1:4*(x)) = sampinex1(x,1:4);
        end
        
        %�����������
        expectout = traind_s( i + inbox,:);



        %���������Լ����ļ���
        
        %�������
        hid_input = sampinex2*net.w1; 
        hid_out = tansig(hid_input); 
        
        %��������
        ou_input1 = [hid_out,1];
        ou_input2 =ou_input1* net.w2 ;
        out_out = tansig(ou_input2);  
        
        %�������
        outrec(:,i) = out_out'; 
        
        %����ÿ�����
        err = expectout - out_out;

        %������
        sse = sumsqr(err);
        %��¼���
        errrec(i,:) = sse;
        
        wuchajilu(i,:) = err;
        %% ���ķ��򴫲�
        if i < 2000
            
            for ii=1:length(hid_out )+1 
                for j=1:4
                    ipu1(j)=err(j);     % �ֲ��ݶ�
                    % �������������֮��ĵ�����
                    delta1(ii,j) = eta.*ipu1(j).*ou_input1(ii); 
                end
            end

                for m=1:40
                    for ii=1:length(hid_out)
                        % �ֲ��ݶ�
                        ipu2(ii)=hid_out(ii).*(1-hid_out(ii)).*sum(ipu1.*net.w2(m,:));
                        % ������������֮��ĵ�����
                        delta2(m,ii)= eta.*ou_input1(m).*ipu2(ii);
                    end
                end
        else
            
            a = wuchajilu(i,:)-wuchajilu(i-1,:);
            b = (net.w2(:,1)-net.w2old(:,1))';
            for zzli = 1:4
                for zzlj = 1:41
                    J(zzli,zzlj) = a(zzli)/b(zzlj);
                end
            end
            
            if wuchajilu(i,1)-wuchajilu(i-1,1)<0
                uk = uk/1.3
            else
                uk = 1.3*uk
                
            end
            
            delta1   =    (     0.05*J'*J  +   1*eye(41)     )  \  J' *  wuchajilu(i,1)   ;
            delta2 = zeros(40);
            

            
        end
        
                net.w2old = net.w2;
                net.w1old = net.w1;

           if i==1 
            net.w2 =  net.w2+eta * delta1;
            net.w1=  net.w1 + eta * delta2;
           else
             net.w2 =  net.w2+eta * delta1 + mc * dWEXold;
             net.w1 = net.w1 + eta * delta2 +  mc * dwexold;
           end

           dWEXold = delta1;
           dwexold = delta2;

    end
    bianshu
    sse = 0;
        wuchajilu(bianshu,:) = err;
        for count = 1:394
             sse = sse+abs(wuchajilu(count,1))+abs(wuchajilu(count,2))+abs(wuchajilu(count,3))+abs(wuchajilu(count,4));
             end
        sse 
         
        if(bianshu == 100)
            i = i;
                zzl = figure(1);
    for zzl = 1:4
        subplot(2,2,zzl);
        plot([11:10+404],outrec(zzl,:),'r');
        hold on 
        plot(traind_s(:,zzl),'b');
        hold off 
    end
        end



    %% ����
    %%figure(1)
    % axis on
    % hold on
    % grid
    % [o,q] = size(errrec);
    %%plot(1:q,errrec,'b-','linewidth',1.5);
    %%axis([0 12952 0 30]);

    %% �˲�
%     jisuanwucha;
%     if (sume<20)
%     break;
%     end 

end

  i= i;
  



