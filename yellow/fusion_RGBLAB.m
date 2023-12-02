clc,close all,clear all;
%颜色修正图和同态滤波图在LAB空间融合
image=imread('E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowbuddhacorrection.jpg');%颜色修正图
Limage=imread('E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowbuddhalvbo.jpg');%对比度增强图
rgb=im2double(image);
figure(1)
imshow(image);title('原图像')
% L=im2double(Limage)*100;
% rgb2lab转换
lab = rgb2lab(rgb,'WhitePoint','d50');
l = lab(:,:,1);%通道分离
a = lab(:,:,2);
b = lab(:,:,3);
figure(2),
imshow(lab(:,:,1),[0 100]);%各通道值显示，因为查出正常显示范围，需要添加范围
title('L通道')
figure(3),
imshow(lab(:,:,2),[-127 128]);title('a通道')
figure(4),
imshow(lab(:,:,3),[-127 128]);title('b通道')


L=im2double(Limage)*100;
%lab通道合并
LAB=cat(3, L,a,b);
cform1=makecform('lab2srgb');   %Lab转到rgb
rgb1 = applycform(LAB, cform1);
 figure(5);
imshow(rgb1);title('灰度L+ab通道+各通道直方图均衡化');
%  imwrite(rgb1,'E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowbuddhadehaze.jpg');

