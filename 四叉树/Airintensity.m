function [count,Airintensity]=Airintensity(img)
%% 返回值 count:偏色通道号（1,2,3）;Airintensity 三通道天空光取值数组
%%
[hei, wid, ~] = size(img);
image = img;
figure, imshow(img); hold on; 
line([1 wid], [hei/2 hei/2], 'Color', 'r','LineWidth',1.5);
line([wid/2 wid/2], [1 hei], 'Color', 'r','LineWidth',1.5);
yoffset = 1;
xoffset = 1;
label = 0.01*hei*wid;  %  迭代停止条件（图像像素大小的1%）
% label = 9;
while(hei * wid > label)
   
    img_ul = img(1: hei/2, 1: wid/2, :); 
    img_ur = img(1:hei/2, wid/2+1:wid, :);
    img_ll = img(hei/2+1:hei, 1:wid/2, :);
    img_lr = img(hei/2+1:hei, wid/2+1:wid, :);
    N = numel(img_ul);

    %公式（5）为寻找条件  
    %目标函数是均值减方差+第二部分（一共三部分）
    dpscore(1) = mean2(img_ul(:,:,1)) - std2(img_ul(:,:,1));
    dpscore(2) = mean2(img_ul(:,:,2)) - std2(img_ul(:,:,2));
    dpscore(3) = mean2(img_ul(:,:,3)) - std2(img_ul(:,:,3));
    afscore(1) = sum(dpscore(:))+1/N.*sum(sum(img_ul(:,:,3)+img_ul(:,:,2)-2.*img_ul(:,:,1)));
    
    dpscore(1) = mean2(img_ur(:,:,1)) - std2(img_ur(:,:,1));
    dpscore(2) = mean2(img_ur(:,:,2)) - std2(img_ur(:,:,2));
    dpscore(3) = mean2(img_ur(:,:,3)) - std2(img_ur(:,:,3));
    afscore(2) = sum(dpscore(:))+1/N.*sum(sum(img_ur(:,:,3)+img_ur(:,:,2)-2.*img_ur(:,:,1)));

    dpscore(1) = mean2(img_ll(:,:,1)) - std2(img_ll(:,:,1));
    dpscore(2) = mean2(img_ll(:,:,2)) - std2(img_ll(:,:,2));
    dpscore(3) = mean2(img_ll(:,:,3)) - std2(img_ll(:,:,3));
    afscore(3) = sum(dpscore(:))+1/N.*sum(sum(img_ll(:,:,3)+img_ll(:,:,2)-2.*img_ll(:,:,1)));
    
    dpscore(1) = mean2(img_lr(:,:,1)) - std2(img_lr(:,:,1));
    dpscore(2) = mean2(img_lr(:,:,2)) - std2(img_lr(:,:,2));
    dpscore(3) = mean2(img_lr(:,:,3)) - std2(img_lr(:,:,3));
    afscore(4) = sum(dpscore(:))+1/N.*sum(sum(img_lr(:,:,3)+img_lr(:,:,2)-2.*img_lr(:,:,1)));
    
    % 找满足条件的背景区域 img为找寻结果
    [~, maxind] = max(afscore);
    clear img;           
    if(maxind == 1)
      img = img_ul;
      [hei, wid, ~] = size(img_ul);
    elseif(maxind == 2)
      img = img_ur;
      [hei, wid, ~] = size(img_ur);
      xoffset = xoffset + wid;
    elseif(maxind == 3)
      img = img_ll;
      [hei, wid, ~] = size(img_ll);
      yoffset = yoffset + hei;
    elseif(maxind == 4)
      img = img_lr;
      [hei, wid, ~] = size(img_lr);
      yoffset = yoffset + hei;
      xoffset = xoffset + wid;
    end
    linehandle1 = line([xoffset, xoffset+wid], [yoffset+hei/2, yoffset+hei/2], 'Color', 'r','LineWidth',1.5);
    linehandle2 = line([xoffset+wid/2, xoffset+wid/2], [yoffset, yoffset+hei], 'Color', 'r','LineWidth',1.5);
end
delete(linehandle1); delete(linehandle2);
%% img为最后被选中的区域 选中区域后坐标取了这个区域的中点
text(xoffset+wid/2,yoffset+hei/2,'o','color','r');
%% 求取img选中区域各个通道的均值并显示
Airintensityr=img(:,:,1);averager = mean(mean(Airintensityr));X = ['选取区域红色通道的均值为 ',num2str(averager)];disp(X);
Airintensityg=img(:,:,2);averageg = mean(mean(Airintensityg));Y = ['选取区域绿色通道的均值为 ',num2str(averageg)];disp(Y);
Airintensityb=img(:,:,3);averageb = mean(mean(Airintensityb));Z = ['选取区域蓝色通道的均值为 ',num2str(averageb)];disp(Z);
 [~, count] = max([averager,averageg,averageb]);
Airintensity = image(xoffset+wid/2,yoffset+hei/2,:);
end