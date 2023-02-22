
z = f_distance;
Nx = 2000;
Ny = 2000;
dm = 0.0024;
Nm = Nx + Ny;
Lm = Nm*dm;
d = 2.5;
r1 = 0.30;
M = (1+d/z);
ri = r1*M;
bi = pi/(r1)^2;
lambda = (440:10:700)*1e-6; 
w = ones(1,length(lambda));
[xm,ym] = meshgrid(-Lm/2:dm:Lm/2-dm);  
mask = 0.5+0.5*square(bi*(xm.^2+ym.^2)+pi/2);  
sup = padarray(ones(Nm-20,Nm-20),[10,10]);
mask = mask.*sup;
f_max = 1 / dm;
du = f_max / Nm;
[u,v] = meshgrid(-Nx*du:du:Ny*du-du);
Im = 0;
for k = 1:length(lambda)
    wave = lambda(k);
    sphere = exp(1i*2*pi*sqrt(xm.^2+ym.^2+z.^2)/wave)/(1i*wave*z); 
    U = mask.*sphere;
    FU = fftshift(fft2(fftshift(U)));
    Ha = exp(1i*2*pi*d*sqrt(1/wave^2-u.^2-v.^2)); 
    U1 = fftshift(ifft2(fftshift(FU.*Ha)));
    Im = w(k)*Im + (U1.*conj(U1));
end
Im = Im(Nm/4:Nm*3/4-1,Nm/4:Nm*3/4-1);
Im = mat2gray(Im);
figure; imagesc(Im);
H = fftshift(fft2(fftshift(real(Im))));

%%

[m,n] = size(I);

HT = conj(H);

% set function handle
ft = @(x) fftshift(fft2(fftshift(x)));
ift = @(x) fftshift(ifft2(fftshift(x)));
C = @(x) conv_H(x,H);
CT = @(x) conv_H(x,HT);

Hh = fftshift(psf2otf([1 -1],[m,n]));HhT = conj(Hh);   
Hv = fftshift(psf2otf([1;-1],[m,n]));HvT = conj(Hv);
Ih = Hh.*HhT;
Iv = Hv.*HvT;
Dh = @(x) conv_H(x,Hh);
Dv = @(x) conv_H(x,Hv);

%% Initialization
x = 0;
bh = Dh(I);
bv = Dv(I);
zh = bh;
zv = bv;

lambda = 0.5;
mu = 1e6;
eta = 0.07;
beta = 0.5;
y1 = 0;
y2 = 0;
y3 = 0;
y4 = 0;
wh = 0;
wv = 0;
denominator = mu*(Ih + Iv) + eta*(Ih + Iv).*abs(H).^2;

N = 20;     % Number of Iterations

figure('OuterPosition',[400 300 800 400]), 
subplot(1,2,1),h = imagesc(x);
t = title('Iteration = 0');
axis image off
colormap gray

subplot(1,2,2),an = animatedline;
xlim([0,N])

tic
%% begin iteration
for k = 1:N
    
    % min x
    numerator = mu*HhT.*ft(wh + y1/mu) + ...
                mu*HvT.*ft(wv + y2/mu) + ...
                eta*HT.*HhT.*ft(zh + y3/eta) + ...
                eta*HT.*HvT.*ft(zv + y4/eta);
    Fx = (numerator + eps)./(denominator + eps);
    x = ift(Fx);
    
    % min w
    sh = ift(Hh.*Fx) - y1/mu;
    sv = ift(Hv.*Fx) - y2/mu;
    [wh,wv] = soft_isotropic(sh,sv,mu);
    
    % min z
    DCxh = ift(Hh.*H.*ft(x));
    DCxv = ift(Hv.*H.*ft(x));
    th = DCxh - y3/eta;
    tv = DCxv - y4/eta;
    zh = (lambda*bh + eta*th)/(lambda + eta);
    zv = (lambda*bv + eta*tv)/(lambda + eta);
    zh = (zh - th) + th;
    zv = (zv - tv) + tv;    
    
    % updata y1,y2,y3,y4
    y1 = y1 - beta*mu*(Dh(x) - wh);
    y2 = y2 - beta*mu*(Dv(x) - wv);
    y3 = y3 - beta*eta*(DCxh - zh);
    y4 = y4 - beta*eta*(DCxv - zv);

    % display
    x_disp = x;
    set(h,'CData',real(x_disp));
    set(t,'string',['Iteration = ',num2str(k)])
    
    resid = abs(Dh(C(x) - I)).^2 + abs(Dv(C(x) - I)).^2;
    resid = sum(resid,'all');
    regul = sqrt(abs(Dh(x)).^2 + abs(Dv(x)).^2);
    regul = sum(regul,'all');
    objective = regul + lambda*resid/2;
    
    addpoints(an,k,objective);
    drawnow
    
end
toc

figure,imshow(x_disp,[]);title('reconstruction')
%% functions
function y = conv_H(x,H)

    Fx = fftshift(fft2(fftshift(x)));
    y = fftshift(ifft2(fftshift(Fx.*H)));

end

function [wh,wv] = soft_isotropic(sh,sv,eta)

    s = sqrt(abs(sh).^2 + abs(sv).^2);
    wh = (sh./s).*max(s - 1/eta,0);
    wv = (sv./s).*max(s - 1/eta,0);

end