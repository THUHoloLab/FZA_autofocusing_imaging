clearvars; clc; close all
addpath('./Function');

%%
img = im2double(imread('./source/experiment/hololab.bmp'));
[Row,Col] = size(img);
img = rot90(img,2);

Xc = Col/2+1;
Yc = Row/2+1;

di = 2.5;             
dp = 0.0024;       
r1 = 0.3;          
Nx = 2000;           
Ny = 2000;

theta = 0.6;

I = img(Yc-Ny/2:Yc+Ny/2-1,Xc-Nx/2:Xc+Nx/2-1);
figure; imagesc(I(250:1700,350:1800)); colormap gray

fu_max = 0.5 / dp;
fv_max = 0.5 / dp;
du = 2*fu_max / Nx;
dv = 2*fv_max / Ny;
[u,v] = meshgrid(-fu_max:du:fu_max-du,-fv_max:dv:fv_max-dv);

%%
start = 100;            %start distance
over = 300;             %over distance
iterations = 100;       %the number of iterations
r_1 = (1+di/start)*r1;
r_2 = (1+di/over)*r1;

iteration = 0;
initial_step = abs(over - start)/iterations;
step = initial_step;
z = start;

wtn = zeros(1,iterations);
gnorm = zeros(1,iterations);
tog = zeros(1,iterations);
x = zeros(1,iterations);

for iteration = 1:iterations
    z = z + step;
    M = di/z;
    r = (1+M)*r1;
    x(iteration) = z;
    Im = set_image(r,I,dp);    
    Im = Gauss_1(Im,300,1);
    gnorm(iteration) = GNORM(Im);
    tog(iteration) = ToG(Im);

end
%%
gnorm = mapminmax(gnorm,0,1);
tog = mapminmax(tog,0,1);
wtn = theta*gnorm + (1-theta)*tog;
figure; plot(x,wtn,'-*','MarkerIndices',1:10:iterations); 
legend('W-T-N');
xlabel('Estimated distance (mm)');
ylabel('Normalized metric');
axis([-inf inf 0 1.1])
box off
[y,x_wtn] = max(wtn);
f_distance = x_wtn*step + start;
M = di/f_distance;
r = (1+M)*r1;
im = set_image(r,I,dp);
figure; imagesc(im); colormap gray; title('W-T-N');




%%
function Im = set_image(r,Im,dp)
    [Ny,Nx] = size(Im);
    ri = r;  
    fu_max = 0.5 / dp;
    fv_max = 0.5 / dp;
    du = 2*fu_max / Nx;
    dv = 2*fv_max / Ny;
    [u,v] = meshgrid(-fu_max:du:fu_max-du,-fv_max:dv:fv_max-dv);
    H = 1i*exp(-1i*(pi*ri^2)*(u.^2 + v.^2));  % fresnel transfer function
    Or = MyAdjointOperatorPropagation(Im,H);
    Im = real(Or);
end




