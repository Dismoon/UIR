clc,close all,clear all;
%��ɫ����ͼ��̬ͬ�˲�ͼ��LAB�ռ��ں�
image=imread('E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\yellowbuddhacorrection.jpg');%��ɫ����ͼ
Limage=imread('E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\yellowbuddhalvbo.jpg');%�Աȶ���ǿͼ
rgb=im2double(image);
figure(1)
imshow(image);title('ԭͼ��')
% L=im2double(Limage)*100;
% rgb2labת��
lab = rgb2lab(rgb,'WhitePoint','d50');
l = lab(:,:,1);%ͨ������
a = lab(:,:,2);
b = lab(:,:,3);
figure(2),
imshow(lab(:,:,1),[0 100]);%��ͨ��ֵ��ʾ����Ϊ���������ʾ��Χ����Ҫ��ӷ�Χ
title('Lͨ��')
figure(3),
imshow(lab(:,:,2),[-127 128]);title('aͨ��')
figure(4),
imshow(lab(:,:,3),[-127 128]);title('bͨ��')


L=im2double(Limage)*100;
%labͨ���ϲ�
LAB=cat(3, L,a,b);
cform1=makecform('lab2srgb');   %Labת��rgb
rgb1 = applycform(LAB, cform1);
 figure(5);
imshow(rgb1);title('�Ҷ�L+abͨ��+��ͨ��ֱ��ͼ���⻯');
%  imwrite(rgb1,'E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\yellowbuddhadehaze.jpg');

