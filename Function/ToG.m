%% ToG
function [out] = ToG(in)
[Nx,Ny]=size(in);
k = 0;
l = gradient(in);
for ii=2:(Nx-1)
    for jj=2:(Ny-1)
        k = k+sqrt((in(ii,jj)-in(ii-1,jj))^2+(in(ii,jj)-in(ii,jj-1))^2);
    end
end
m = std(l,0,"all");
out = sqrt(m*Nx*Ny/k);
end