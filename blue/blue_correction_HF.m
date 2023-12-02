clc;
clear all;
close all;
image=imread('E:\F盘\图像处理\水下图像\Reg_ Research Help\bluebuddha.jpg');
%批量处理

figure();
imshow(image);
gray=rgb2gray(image);
Gray1=im2double(gray);
%同态滤波
rH =1.2;
rL = 0.5;
c = 0.6;%介于rH和rL之间
D0 =2000;
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
img_after1=adapthisteq(img_after1);
% img_after1=1.2*(img_after1);
% imwrite(img_after1,'E:\F盘\图像处理\自编\灰度滤波+通道补偿\img_after1.jpg');
% figure()
% subplot(1,2,1), imshow(Gray1), title('灰度图');
% subplot(1,2,2), imshow(img_after1), title('滤波图');
%  imwrite(img_after1,'C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\blue2\160lvbo.jpg');
% % %颜色校正

%%R通道
R1 = image(:,:,1);
R=im2double(R1);
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
% subplot(2,2,3);
% imshow(B)
%  title('B通道');
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
 
% %红色通道补偿
% Irc=R+2*(b-a)*(1-R).*G;%绿通道补偿红通道

Irc=R+1.2*(c-a)*(1-R).*B;%蓝通道补偿红通道
Irc1 = im2uint8(Irc);%有可能值超过1，所以转换到0-255，，uint8(round(I*255)); 

Colorbuchang(:,:,1)=Irc1(:,:,1);
Colorbuchang(:,:,2)=G1  (:,:,1);
Colorbuchang(:,:,3)=B1(:,:,1);

% % %红色通道动态调整
%%R通道
R1 = Colorbuchang(:,:,1);
R=im2double(R1);
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
G1 = Colorbuchang(:,:,2);
G=im2double(G1);

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
B1 = Colorbuchang(:,:,3);
B=im2double(B1);%数据从0~255映射到0~1
% subplot(2,2,3);
% imshow(B)
%  title('B通道');
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
 
a1=8000;
A=(m*n)/a1;
[Ir,Pr]=imhist(R1);%
% %寻找最大值
% [data_max,index]=max(R1(:));

%
Ir(Ir<A) = 0;
[Lrcmin,~] = find(Ir>0,1,'first');%参考基于色彩平衡与融合的水下图像增强技术（哈工程）
[Lrcmax,~] = find(Ir>0,1,'last');
Drc=Lrcmax-Lrcmin;
%增益
Kr=255/Drc;
%通道拉伸
Irrc=Kr*(R1-Lrcmin);

%绿色
[Ig,Pg]=imhist(G1);%
% %寻找最大值
% [data_max,index]=max(R1(:));
%
Ig(Ig<A) = 0;
[Lgcmin,~] = find(Ig>0,1,'first');
[Lgcmax,~] = find(Ig>0,1,'last');
Dgc=Lgcmax-Lgcmin;
%增益
Kg=255/Dgc;
%通道拉伸
Iggc=Kg*(G1-Lgcmin);
% figure()
% imshow(Iggc);
% figure()
% imshow(G1);
%蓝色
[Ib,Pb]=imhist(B1);%
% %寻找最大值
% [data_max,index]=max(R1(:));
%
Ib(Ib<A) = 0;

[Lbcmin,~] = find(Ib>0,1,'first');
[Lbcmax,~] = find(Ib>0,1,'last');
Dbc=Lbcmax-Lbcmin;
%增益
Kb=255/Dbc;
%通道拉伸
Ibbc=Kb*(B1-Lbcmin);

% % 三通道合成彩色图像
Colorxiu(:,:,1)=Irrc(:,:,1);
Colorxiu(:,:,2)=Iggc(:,:,1);
Colorxiu(:,:,3)=Ibbc(:,:,1);

rgb=im2double(Colorxiu);
% figure()
% imshow(Colorxiu);title('颜色修正图像')
%
% L=im2double(Limage)*100;
% rgb2lab转换
lab = rgb2lab(rgb,'WhitePoint','d50');
l = lab(:,:,1);%通道分离
a_lab = lab(:,:,2);
b_lab = lab(:,:,3);

L=im2double(img_after1)*100;
%lab通道合并
LAB=cat(3, 1.2*L,a_lab,b_lab);
cform1=makecform('lab2srgb');   %Lab转到rgb
rgb1 = applycform(LAB, cform1);
 figure();
imshow(rgb1);title('灰度L+ab通道+各通道直方图均衡化');
%  imwrite(rgb1,'C:\Users\Administrator\Desktop\水下彩色图像校正数据（专利）\blue2\RUIE\blue_92_our1.jpg');
