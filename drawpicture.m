clc 
clear all 

load out2.mat
load train.mat
outrec = outrec';
ml = [2.35500000000000,0.110000000000000,6.50000000000000,323.500000000000];
mm = [10.9945239746905,69.0063006300631,-6.66869763899466,-0.319387271030040];
for j = 1:4
    for i = 1:404
        outrec(i,j) = outrec(i,j)*ml(j)+mm(j);

    end
end

for j = 1:4
    for i = 1:404
        traind_s(i,j) = traind_s(i,j)*ml(j)+mm(j);

    end
end


xunlian = 250;
yuce = 152;

a = 1:xunlian;
b = xunlian+1:xunlian+yuce;

figure(1)
%% ÄâºÏÇúÏß

title('fitted');

zzl = 1;
subplot(2,2,zzl);
plot(a,traind_s(1:xunlian,zzl),'b');
hold on 
plot(a+10,outrec(1:xunlian,zzl),'r');
hold off 
xlabel('time');
ylabel('temperature');
legend('measured value','fitted value');

        
zzl =2;
subplot(2,2,zzl);
plot(a,traind_s(1:xunlian,zzl),'b');
hold on 
plot(a+10,outrec(1:xunlian,zzl),'r');
hold off 
xlabel('time');
ylabel('pH');
legend('measured value','fitted value');

zzl = 3;
subplot(2,2,zzl);
plot(a,traind_s(1:xunlian,zzl),'b');
hold on 
plot(a+10,outrec(1:xunlian,zzl),'r');
hold off 
xlabel('time');
ylabel('Do');
legend('measured value','fitted value');

zzl = 4;
subplot(2,2,zzl);
plot(a,traind_s(1:xunlian,zzl),'b');
hold on 
plot(a+10,outrec(1:xunlian,zzl),'r');
hold off 
xlabel('time');
ylabel('ORP');
legend('measured value','fitted value');

%% Ô¤²âÇúÏß
figure(2)

zzl = 1;
subplot(2,2,zzl);

plot(b,traind_s(xunlian+1:xunlian+yuce,zzl),'b');
hold on 
plot(b+10,outrec(xunlian+1:xunlian+yuce,zzl),'r');
hold off 

xlabel('time');
ylabel('temperature');
legend('measured value','predicted value');

zzl = 2;
subplot(2,2,zzl);
plot(b,traind_s(xunlian+1:xunlian+yuce,zzl),'b');
hold on 
plot(b+10,outrec(xunlian+1:xunlian+yuce,zzl),'r');
hold off 
xlabel('time');
ylabel('pH');
legend('measured value','predicted value');

zzl = 3;
subplot(2,2,zzl);
plot(b,traind_s(xunlian+1:xunlian+yuce,zzl),'b');
hold on 
plot(b+10,outrec(xunlian+1:xunlian+yuce,zzl),'r');
hold off 
xlabel('time');
ylabel('DO');
legend('measured value','predicted value');

zzl = 4;
subplot(2,2,zzl);
plot(b,traind_s(xunlian+1:xunlian+yuce,zzl),'b');
hold on 
plot(b+10,outrec(xunlian+1:xunlian+yuce,zzl),'r');
hold off 
xlabel('time');
ylabel('ORP');
legend('measured value','predicted value');