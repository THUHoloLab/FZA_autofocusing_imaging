%% 方差
function [out] = VAR(in)
[Nx,Ny]=size(in);
k=0;
in  = Gauss_1(in,200,1);
m = mean(mean(in));
for ii=2:(Nx-1)
    for jj=2:(Ny-1)
        k=k+(in(ii,jj) - m)^2;
    end
end
out = k/(Nx*Ny);
end