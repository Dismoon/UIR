clc;
clear all;
close all;
%��ɫͼ��ת���Ҷ�ͼ������̬ͬ�˲�
image=imread('E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\yellowcolorchart.jpg');
figure();
imshow(image);
title('ԭʼͼ��');
gray=rgb2gray(image);
figure,imshow(gray);
Gray1=im2double(gray);

%̬ͬ�˲�
rH =1.1;
rL = 0.5;
c = 0.6;%����rH��rL֮��
D0 =2000;%2000

[M, N] = size(Gray1);

%ȡ����
img_log1 = log(double(Gray1) + 1);


%ƽ�Ƶ����ģ��ж�������ָ������
img_py1 = zeros(M, N);
for i = 1:M
   for j= 1:N
       if mod(i+j, 2) == 0
           img_py1(i,j) = img_log1(i, j);
       else
           img_py1(i,j) = -1 * img_log1(i, j);
       end
   end
end

% �������ͼ����и���Ҷ�任
img_py_fft1 = fft2(img_py1);

%̬ͬ�˲�����
img_tt1 = zeros(M, N);


u=floor(M/2);%���ĵ�����
v=floor(N/2);  

for i = 1:M
   for j =1:N
        D(i,j) =sqrt( ((i-u).^2+(j-v).^2));
        img_tt1(i, j) =(rH-rL) .* (1-(exp(c*(-D(i,j)^2./D0^2)))) + rL;
%         img_tt(i, j)=(rH-rL).*(1/(1+(D0/c*(D(i,j).*1)^2)))+rL;%������˹�˲���
%             img_tt(i, j)=(rH-rL)*(1/(1+(c*D0/(D(i,j))).^2))+rL;%������˹�˲���
%             img_tt(i,j)=(rH-rL) .* (exp((-D0^2)./D(i,j).^2))+ rL ;%ָ���ʹ��ݺ���
            
   end
end

% figure()
% mesh(D);%�˲���ʾ��ͼ


%�˲�
img_temp1 =   img_py_fft1.*img_tt1;

%���任,ȡʵ��������ֵ
img_temp1 = abs(real(ifft2(img_temp1)));

%ָ����
img_temp1 = exp(img_temp1) - 1;


%��һ������
max_num1 = max(img_temp1(:));
min_num1 = min(img_temp1(:));
range1 = max_num1 - min_num1;
img_after1 = zeros(M,N,'uint8');
for i = 1 : M
    for j = 1 : N
        img_after1(i,j) = uint8(255 * (img_temp1(i, j)-min_num1) / range1);
    end
end
%ԭͼ���˲���ͼ��Ա�
img_after2=adapthisteq(img_after1);
%  imwrite(img_after2,'E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\yellowcolorchartlvbo.jpg');
figure();
subplot(1,3,1), imshow(Gray1), title('�Ҷ�ͼ');
subplot(1,3,2), imshow(img_after1), title('�˲�ͼ');
subplot(1,3,3), imshow(img_after2), title('clahe�˲�ͼ');