function output = color_enhancement(image)
    % 颜色增强函数
    % 输入参数：
    % image - 输入图像
    
    % 颜色通道补偿
    output = color_channel_compensation(image);
end

function output = color_channel_compensation(image, win_size, k, landa)
    % 颜色通道补偿函数
    % 输入参数：
    % image - 输入图像
    % win_size - 窗口大小（默认20）
    % k - 补偿系数（默认1）
    % landa - 补偿系数（默认0.5）
    
    if ~exist('win_size','var')
        win_size = 20;
    end
    
    if ~exist('k','var')
        k = 1;
    end
    
    if ~exist('landa','var')
        landa = 0.5;
    end
    
    [m,n,~] = size(image);
    im_RGB = im2double(image);
    
    % 转换到LAB颜色空间
    cform = makecform('srgb2lab');
    lab = applycform(image,cform);
    lab = lab2double(lab);
    I_L = lab(:,:,1);
    I_a = lab(:,:,2);
    I_b = lab(:,:,3);
    
    % 高斯滤波
    gausFilter = fspecial('gaussian',[20,20],1);
    gaus_a = imfilter(I_a,gausFilter,'replicate');
    gaus_b = imfilter(I_b,gausFilter,'replicate');
    
    % 图像扩展
    Im = zeros(m+(win_size-1),n+(win_size-1),3);
    radius_size = floor(win_size/2);
    Im(radius_size:(m+radius_size-1),radius_size:(n+radius_size-1),:) = im_RGB;
    
    % 使用矩阵运算计算局部均值
    im_mean = zeros(m,n);
    for i = 1:3
        temp = conv2(Im(:,:,i), ones(win_size,win_size)/(win_size*win_size), 'valid');
        im_mean = im_mean + temp;
    end
    im_mean = im_mean/3;
    
    % 使用矩阵运算进行颜色补偿
    M = ones(m,n);
    M(im_mean > 0.85) = 0;
    
    Ic_a = I_a - k*M.*gaus_a;
    Ic_b = I_b - landa*M.*gaus_b;
    
    % 定义形态学操作核
    kernel1 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
               1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    
    kernel2 = [1,1,1,1,1
               1,0,0,0,1
               1,0,0,0,1
               1,0,0,0,1
               1,1,1,1,1];
    
    % 应用形态学变换
    lab(:,:,1) = MTHT(I_L,10,kernel1,kernel2,0.4);
    lab(:,:,2) = Ic_a;
    lab(:,:,3) = Ic_b;
    
    % 转换回RGB颜色空间
    output = lab2rgb(lab);
end

function output = MTHT(image, time, block1, block2, weight)
    % 多尺度形态学变换函数（辅助函数）
    if ~exist('time','var')
        time = 10;
    end
    
    if ~exist('block1','var')
        block1 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,0,0,0,0,0,0,0,0,0,0,0,0,0,1
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    end
    
    if ~exist('block2','var')
        block2 = [1,1,1,1,1
                 1,0,0,0,1
                 1,0,0,0,1
                 1,0,0,0,1
                 1,1,1,1,1];
    end
    
    if ~exist('weight','var')
        weight = 0.1;
    end
    
    [m,n,d] = size(image);
    
    THNM = zeros(m,n,d);
    BHNM = zeros(m,n,d);
    THNSM = zeros(m,n,d);
    BHNSM = zeros(m,n,d);
    
    for i = 1:time
        THN = image-imclose((imopen(image,block1)),block2);
        BHN = imopen((imclose(image,block1)),block2)-image;
        
        THNM = THNM+THN;
        BHNM = BHNM+BHN;
        
        if i > 1
            THNS = THN-T;
            BHNS = BHN-B;
            
            THNSM = THNSM+THNS;
            BHNSM = BHNSM+BHNS;
            
            T = THN;
            B = BHN;
        else
            T = THN;
            B = BHN;
        end
    end
    
    if i == 1
        detail = THN-BHN;
        output = image+detail;
    else
        THNM = weight*THNM;
        BHNM = weight*BHNM;
        THNSM = weight*THNSM;
        BHNSM = weight*BHNSM;
        
        detail = (THNM+THNSM)-(BHNM+BHNSM);
        output = image+detail;
    end
end 