clc;
clear all;
close all;
%%======利用解析法进行直接反投影======%%
N = 2048;
I = phantom(N);
delta = pi/180;%角度增量
theta = 0:1:179;%投影角度
theta_num = length(theta);%总共有180个角度
%%=====生产投影数据=====%%
P = radon(I,theta);
[mm,nn] = size(P);%计算投影矩阵的行和列
e = floor((mm-N-1)/2+1)+1;%投影数据的默认投影中心为floor(size(I)+1/2)
P = P(e:N+e-1,:);%截取中心N点数据,因投影数据较多，含无用数据
P1=reshape(P,N,theta_num);
%%=====反投影重建=====%%
rec = medfuncBackprojection(theta_num,N,P1,delta);
%%=====投影结果显示=====%%
figure;
subplot(1,2,1),imshow(I,[]),title('原始图像');
subplot(1,2,2),imshow(rec,[]),title('重建后的图像');
%%=====子程序=====%%
function rec = medfuncBackprojection(theta_num,N,R1,delta)
%Backprojection reconstruction function
%---------------
%输入参数：
%theta_num：投影角度个数
%N：图像大小，探测器通道数
%R1：投影数据矩阵（N*theta_num）
%delta:角度增量
%：反投影重建图像矩阵
%==============================================================%
rec = zeros(N);%用于存储重建后的像素值
for m = 1:theta_num
    pm=R1(:,m);%取某一角度的投影数据,R1投影数据矩阵
    cm = (N/2)*(1-cos((m-1)*delta)-sin((m-1)*delta));
    for k1=1:N
        for k2 = 1:N
            %以下是射束计算，需要注意的是射数编号n取值范围为1~N-1
            xrm= cm+(k2-1)*cos((m-1)*delta)+(k1-1)*sin((m-1)*delta);
            n=floor(xrm);%射束编号（整数部分）
            t=xrm-n;%小数部分
            n=max(1,n);n=min(n,N-1);%限定n范围为1~N-1
            p=(1-t)*pm(n)+t*pm(n+1);%线性内插
            %rec(k1,k2)=rec(k1,k2)+p;%反投影，图像需要翻转180°
            rec(N+1-k1,k2)=rec(N+1-k1,k2)+p;%反投影，图像需要翻转90°
        end
    end       
end
end