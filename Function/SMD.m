%% SMD（灰度方差）函数
function [out] = SMD(in)
[Nx,Ny]=size(in);
k=0; 
in_1 = mat2gray(in);
for ii=2:(Nx-1)
    for jj=2:(Ny-1)
        k=k+(abs(in_1(ii,jj)-in(ii,jj-1))+abs(in(ii,jj)-in(ii+1,jj)));
    end
end
out = k;
end