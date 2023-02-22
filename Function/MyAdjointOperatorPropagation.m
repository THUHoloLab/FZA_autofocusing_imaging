function Or = MyAdjointOperatorPropagation(I,H)

FI = fftshift(fft2(fftshift(I)));
%fftshift 将零频分量移到频谱中心
% Or = fftshift(ifft2(fftshift(FI.*conj(H))));
Or = fftshift(ifft2(fftshift(FI./H)));
%ifft2 二维快速傅里叶逆变换

Or = real(Or);  %取实部
% Or = Or/max(Or(:));

end