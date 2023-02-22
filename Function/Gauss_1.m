function FFT = Gauss_1(im,sigma,model)

[N,M] = size(im);
tic
if sigma < 0
    sigma = 100;
end
F = fftshift(fft2(fftshift(im)));
x=-M/2:M/2-1;
y=-N/2:N/2-1;
u1 = 0;          %均值
u2 = 0;
l = N/M;
sigma1 = sigma;    %方差
sigma2 = sigma*l;
rou = 0;     %相关系数
[X,Y]=meshgrid(x,y); % 产生网格数据并处理
p = 1/(2*pi*sigma1*sigma2*sqrt(1-rou*rou)).*exp(-1/(2*(1-rou^2)).*((X-u1).*(X-u1)/(sigma1*sigma1)-2*rou*(X-u1).*(Y-u2)/(sigma1*sigma2)+(Y-u2).*(Y-u2)/(sigma2*sigma2)));

if model == 1 
    p = p/p(ceil(N/2),ceil(M/2));
end
if model == 2
    p = (1-p/p(N/2,M/2));
end

blur_fft = F.*p.*p.*p.*p;
blur_img = fftshift(ifft2(fftshift(blur_fft)));
FFT = abs(blur_img);
end