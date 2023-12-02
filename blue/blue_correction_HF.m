clc;
clear all;
close all;
image=imread('E:\F��\ͼ����\ˮ��ͼ��\Reg_ Research Help\bluebuddha.jpg');
%��������

figure();
imshow(image);
gray=rgb2gray(image);
Gray1=im2double(gray);
%̬ͬ�˲�
rH =1.2;
rL = 0.5;
c = 0.6;%����rH��rL֮��
D0 =2000;
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
img_after1=adapthisteq(img_after1);
% img_after1=1.2*(img_after1);
% imwrite(img_after1,'E:\F��\ͼ����\�Ա�\�Ҷ��˲�+ͨ������\img_after1.jpg');
% figure()
% subplot(1,2,1), imshow(Gray1), title('�Ҷ�ͼ');
% subplot(1,2,2), imshow(img_after1), title('�˲�ͼ');
%  imwrite(img_after1,'C:\Users\Administrator\Desktop\ˮ�²�ɫͼ��У�����ݣ�ר����\blue2\160lvbo.jpg');
% % %��ɫУ��

%%Rͨ��
R1 = image(:,:,1);
R=im2double(R1);
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
% subplot(2,2,3);
% imshow(B)
%  title('Bͨ��');
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
 
% %��ɫͨ������
% Irc=R+2*(b-a)*(1-R).*G;%��ͨ��������ͨ��

Irc=R+1.2*(c-a)*(1-R).*B;%��ͨ��������ͨ��
Irc1 = im2uint8(Irc);%�п���ֵ����1������ת����0-255����uint8(round(I*255)); 

Colorbuchang(:,:,1)=Irc1(:,:,1);
Colorbuchang(:,:,2)=G1  (:,:,1);
Colorbuchang(:,:,3)=B1(:,:,1);

% % %��ɫͨ����̬����
%%Rͨ��
R1 = Colorbuchang(:,:,1);
R=im2double(R1);
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
G1 = Colorbuchang(:,:,2);
G=im2double(G1);

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
B1 = Colorbuchang(:,:,3);
B=im2double(B1);%���ݴ�0~255ӳ�䵽0~1
% subplot(2,2,3);
% imshow(B)
%  title('Bͨ��');
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
 
a1=8000;
A=(m*n)/a1;
[Ir,Pr]=imhist(R1);%
% %Ѱ�����ֵ
% [data_max,index]=max(R1(:));

%
Ir(Ir<A) = 0;
[Lrcmin,~] = find(Ir>0,1,'first');%�ο�����ɫ��ƽ�����ںϵ�ˮ��ͼ����ǿ�����������̣�
[Lrcmax,~] = find(Ir>0,1,'last');
Drc=Lrcmax-Lrcmin;
%����
Kr=255/Drc;
%ͨ������
Irrc=Kr*(R1-Lrcmin);

%��ɫ
[Ig,Pg]=imhist(G1);%
% %Ѱ�����ֵ
% [data_max,index]=max(R1(:));
%
Ig(Ig<A) = 0;
[Lgcmin,~] = find(Ig>0,1,'first');
[Lgcmax,~] = find(Ig>0,1,'last');
Dgc=Lgcmax-Lgcmin;
%����
Kg=255/Dgc;
%ͨ������
Iggc=Kg*(G1-Lgcmin);
% figure()
% imshow(Iggc);
% figure()
% imshow(G1);
%��ɫ
[Ib,Pb]=imhist(B1);%
% %Ѱ�����ֵ
% [data_max,index]=max(R1(:));
%
Ib(Ib<A) = 0;

[Lbcmin,~] = find(Ib>0,1,'first');
[Lbcmax,~] = find(Ib>0,1,'last');
Dbc=Lbcmax-Lbcmin;
%����
Kb=255/Dbc;
%ͨ������
Ibbc=Kb*(B1-Lbcmin);

% % ��ͨ���ϳɲ�ɫͼ��
Colorxiu(:,:,1)=Irrc(:,:,1);
Colorxiu(:,:,2)=Iggc(:,:,1);
Colorxiu(:,:,3)=Ibbc(:,:,1);

rgb=im2double(Colorxiu);
% figure()
% imshow(Colorxiu);title('��ɫ����ͼ��')
%
% L=im2double(Limage)*100;
% rgb2labת��
lab = rgb2lab(rgb,'WhitePoint','d50');
l = lab(:,:,1);%ͨ������
a_lab = lab(:,:,2);
b_lab = lab(:,:,3);

L=im2double(img_after1)*100;
%labͨ���ϲ�
LAB=cat(3, 1.2*L,a_lab,b_lab);
cform1=makecform('lab2srgb');   %Labת��rgb
rgb1 = applycform(LAB, cform1);
 figure();
imshow(rgb1);title('�Ҷ�L+abͨ��+��ͨ��ֱ��ͼ���⻯');
%  imwrite(rgb1,'C:\Users\Administrator\Desktop\ˮ�²�ɫͼ��У�����ݣ�ר����\blue2\RUIE\blue_92_our1.jpg');
