%% Laplacian 梯度函数
function [out] = LAP(in) 
[Nx,Ny]=size(in);
%in  = Gauss_1(in,180,1);
k=0;
for ii=2:(Nx-1)
    for jj=2:(Ny-1)
        k=k+(in(ii+1,jj)+in(ii-1,jj)+in(ii,jj+1)+in(ii,jj-1)-4*in(ii,jj))^2;
    end
end
out=k;
end