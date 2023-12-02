clc;
clear all;
close all;
%彩色图像转到灰度图，进行同态滤波
image=imread('E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowcolorchart.jpg');
figure();
imshow(image);
title('原始图像');
gray=rgb2gray(image);
figure,imshow(gray);
Gray1=im2double(gray);

%同态滤波
rH =1.1;
rL = 0.5;
c = 0.6;%介于rH和rL之间
D0 =2000;%2000

[M, N] = size(Gray1);

%取对数
img_log1 = log(double(Gray1) + 1);


%平移到中心，判断语句代替指数计算
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

% 对填充后的图像进行傅里叶变换
img_py_fft1 = fft2(img_py1);

%同态滤波函数
img_tt1 = zeros(M, N);


u=floor(M/2);%中心点坐标
v=floor(N/2);  

for i = 1:M
   for j =1:N
        D(i,j) =sqrt( ((i-u).^2+(j-v).^2));
        img_tt1(i, j) =(rH-rL) .* (1-(exp(c*(-D(i,j)^2./D0^2)))) + rL;
%         img_tt(i, j)=(rH-rL).*(1/(1+(D0/c*(D(i,j).*1)^2)))+rL;%巴特沃斯滤波器
%             img_tt(i, j)=(rH-rL)*(1/(1+(c*D0/(D(i,j))).^2))+rL;%巴特沃斯滤波器
%             img_tt(i,j)=(rH-rL) .* (exp((-D0^2)./D(i,j).^2))+ rL ;%指数型传递函数
            
   end
end

% figure()
% mesh(D);%滤波器示意图


%滤波
img_temp1 =   img_py_fft1.*img_tt1;

%反变换,取实部，绝对值
img_temp1 = abs(real(ifft2(img_temp1)));

%指数化
img_temp1 = exp(img_temp1) - 1;


%归一化处理
max_num1 = max(img_temp1(:));
min_num1 = min(img_temp1(:));
range1 = max_num1 - min_num1;
img_after1 = zeros(M,N,'uint8');
for i = 1 : M
    for j = 1 : N
        img_after1(i,j) = uint8(255 * (img_temp1(i, j)-min_num1) / range1);
    end
end
%原图和滤波后图像对比
img_after2=adapthisteq(img_after1);
%  imwrite(img_after2,'E:\F盘\图像处理\水下图像\Reg_ Research Help\yellowcolorchartlvbo.jpg');
figure();
subplot(1,3,1), imshow(Gray1), title('灰度图');
subplot(1,3,2), imshow(img_after1), title('滤波图');
subplot(1,3,3), imshow(img_after2), title('clahe滤波图');