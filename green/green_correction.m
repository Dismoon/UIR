clc;clear;close all;
%绿色偏色处理
% t=imread('E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowcolorchart.jpg'); 
t=imread('C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\green\906_img_.png');

I=t(:,:,1);

[height,width] = size(I);  

subplot(121);

imshow(t),title('原始图像')%显示原始图像   

%对R通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  

%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  

    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  

for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  

end  

%计算累计直方图分布  

c = zeros(1,256);  

c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  

end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  

%对图像进行均衡化

for i = 1:height  

    for j = 1: width  

        Ir(i,j) = c(I(i,j)+1);  

    end  

end  

dis(:,:,1)=Ir;

%subplot(122)  

%imshow(Ir)%显示均衡化后的图像

%对G通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

I=t(:,:,2);

[height,width] = size(I);  

 

%下面使用直方图均衡化进行处理

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  

%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  

    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  

for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  

end  

%计算累计直方图分布  

c = zeros(1,256);  

c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  

end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  

%对图像进行均衡化

for i = 1:height  

    for j = 1: width  

        Ig(i,j) = c(I(i,j)+1);  

    end  

end  

%subplot(122)  

%imshow(Ig)%显示均衡化后的图像

dis(:,:,2)=Ig;

 

%对B通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

I=t(:,:,3);

[height,width] = size(I);  

 

%下面使用直方图均衡化进行处理

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  

%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  

    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  

for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  

end  

%计算累计直方图分布  

c = zeros(1,256);  

c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  

end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  

%对图像进行均衡化

for i = 1:height  

    for j = 1: width  

        Ib(i,j) = c(I(i,j)+1);  

    end  

end  

dis(:,:,3)=Ib; 

%subplot(122)  

%imshow(Ib)%显示均衡化后的图像 

subplot(122);
% imwrite(dis,'E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowcolorchartcorrection.jpg');
imshow(dis),title('各通道直方图均衡化')%显示均衡化后的图像
