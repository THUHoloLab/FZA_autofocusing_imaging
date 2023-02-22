%% Norm 
function [out] = GNORM(in)

[gx,gy] = gradient(in);
Gr = (gx)^2 + (gy)^2;
out = norm(Gr,'fro')/size(Gr,1);

end