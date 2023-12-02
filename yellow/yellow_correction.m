% 黄色图像颜色修正，多蓝色和绿色通道补偿
clc,close all,clear all;
image=imread('C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\yellow\000.jpg');
% [image,rect]=imcrop(I2);%裁剪
% imwrite(image,'C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\blue\112jiequ.jpg');
%%R通道
R1 = image(:,:,1);
R=im2double(R1);

subplot(2,2,1);
imshow(R)
 title('R通道');
%求红色通道灰度平均值
 [m,n]=size(R);
s=0;
for x=1:m
    for y=1:n
        s=s+R(x,y); %求像素值总和 s
    end
end
%所有像素均值
a=mean(mean(R));%先计算列向量均值，再求总均值。
 
%%G通道
G1 = image(:,:,2);
G=im2double(G1);
subplot(2,2,2);
imshow(G);
 title('G通道');
 
 %求绿色通道灰度平均值
 [m,n]=size(G);
k=0;
for x=1:m
    for y=1:n
        k=k+G(x,y); %求像素值总和 s
    end
end
%所有像素均值
b=mean(mean(G));%先计算列向量均值，再求总均值。

%%B通道
B1 = image(:,:,3);
B=im2double(B1);%数据从0~255映射到0~1
subplot(2,2,3);
imshow(B)
 title('B通道');
 %求蓝色通道灰度平均值a
 [m,n]=size(B);
w=0;
for x=1:m
    for y=1:n
        w=w+B(x,y); %求像素值总和 s
    end
end
%所有像素均值
c=mean(mean(B));%先计算列向量均值，再求总均值。
subplot(2,2,4);
imshow(image)
title('原图像');
%绿色通道作为参考通道（在应用中还有些问题。实际为哪个通道值大，哪个作为参考通道）
%红色通道补偿
Irc=R+1.5*(b-a)*(1-R).*G;%绿通道补偿红通道
Irc1 = im2uint8(Irc);%有可能值超过1，所以转换到0-255，，uint8(round(I*255)); 

% % %蓝色通道补偿
Ibc=B+2*(b-c)*(1-B).*G;
Ibc1 = im2uint8(Ibc);%有可能值超过1，所以转换到0-255，，uint8(round(I*255)); 
%图像合成
RGB3(:,:,1)=Irc1 (:,:,1);
RGB3(:,:,2)=G1 (:,:,1);
RGB3(:,:,3)=Ibc1 (:,:,1);
figure();
imshow(RGB3);
title('补偿三通道0-255合成图像');
% imwrite(RGB3,'C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\yellow\748颜色补偿.jpg');
figure()
subplot(121),imshow(image);title('原始图像')
subplot(122),imshow(RGB3);title('通道补偿后图像')

