% ��ɫͼ����ɫ����������ɫ����ɫͨ������
clc,close all,clear all;
image=imread('C:\Users\Administrator\Desktop\ˮ�²�ɫͼ��У�����ݣ�ר����\yellow\000.jpg');
% [image,rect]=imcrop(I2);%�ü�
% imwrite(image,'C:\Users\Administrator\Desktop\ˮ�²�ɫͼ��У�����ݣ�ר����\blue\112jiequ.jpg');
%%Rͨ��
R1 = image(:,:,1);
R=im2double(R1);

subplot(2,2,1);
imshow(R)
 title('Rͨ��');
%���ɫͨ���Ҷ�ƽ��ֵ
 [m,n]=size(R);
s=0;
for x=1:m
    for y=1:n
        s=s+R(x,y); %������ֵ�ܺ� s
    end
end
%�������ؾ�ֵ
a=mean(mean(R));%�ȼ�����������ֵ�������ܾ�ֵ��
 
%%Gͨ��
G1 = image(:,:,2);
G=im2double(G1);
subplot(2,2,2);
imshow(G);
 title('Gͨ��');
 
 %����ɫͨ���Ҷ�ƽ��ֵ
 [m,n]=size(G);
k=0;
for x=1:m
    for y=1:n
        k=k+G(x,y); %������ֵ�ܺ� s
    end
end
%�������ؾ�ֵ
b=mean(mean(G));%�ȼ�����������ֵ�������ܾ�ֵ��

%%Bͨ��
B1 = image(:,:,3);
B=im2double(B1);%���ݴ�0~255ӳ�䵽0~1
subplot(2,2,3);
imshow(B)
 title('Bͨ��');
 %����ɫͨ���Ҷ�ƽ��ֵa
 [m,n]=size(B);
w=0;
for x=1:m
    for y=1:n
        w=w+B(x,y); %������ֵ�ܺ� s
    end
end
%�������ؾ�ֵ
c=mean(mean(B));%�ȼ�����������ֵ�������ܾ�ֵ��
subplot(2,2,4);
imshow(image)
title('ԭͼ��');
%��ɫͨ����Ϊ�ο�ͨ������Ӧ���л���Щ���⡣ʵ��Ϊ�ĸ�ͨ��ֵ���ĸ���Ϊ�ο�ͨ����
%��ɫͨ������
Irc=R+1.5*(b-a)*(1-R).*G;%��ͨ��������ͨ��
Irc1 = im2uint8(Irc);%�п���ֵ����1������ת����0-255����uint8(round(I*255)); 

% % %��ɫͨ������
Ibc=B+2*(b-c)*(1-B).*G;
Ibc1 = im2uint8(Ibc);%�п���ֵ����1������ת����0-255����uint8(round(I*255)); 
%ͼ��ϳ�
RGB3(:,:,1)=Irc1 (:,:,1);
RGB3(:,:,2)=G1 (:,:,1);
RGB3(:,:,3)=Ibc1 (:,:,1);
figure();
imshow(RGB3);
title('������ͨ��0-255�ϳ�ͼ��');
% imwrite(RGB3,'C:\Users\Administrator\Desktop\ˮ�²�ɫͼ��У�����ݣ�ר����\yellow\748��ɫ����.jpg');
figure()
subplot(121),imshow(image);title('ԭʼͼ��')
subplot(122),imshow(RGB3);title('ͨ��������ͼ��')

