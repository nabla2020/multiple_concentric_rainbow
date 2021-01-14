clear; clc;
data=load('wl405nm.txt');  % load the output file
lam0=0.405; % light wavelength, unit: um
nb=1;              %background index is nb=1
R0=5;
d=-R0:R0/10000:-2/3*R0;
ray_index=find(data(:,3)<0);   %find the TIR rays
qx=data(ray_index,2); 
qy=data(ray_index,3);
kx=data(ray_index,4);
ky=data(ray_index,5);
I=data(ray_index,6);
r1=data(ray_index,7);
phase=data(ray_index,8);
d=d(ray_index);
[theta,~]=cart2pol(kx,ky); % to calcuate the ouput angle of light ray
[~,pos_p]=findpeaks(theta);
pos_m=pos_p-1;             
pos_m=[pos_m;length(qy)];
pos_p=[1;pos_p];
for num=1:length(pos_p)
    phase(pos_p(num):pos_m(num))=unwrap(phase(pos_p(num):pos_m(num)));
end
scan=linspace(-pi+0.01,-0.01,10000);   % the range of angle
num=0;
for angle0=scan
    num=num+1;
    for order=1:2;  % the order of the rays
        range_order=pos_p(end+1-order):pos_m(end+1-order);
        I_a(order)=interp1(theta(range_order),I(range_order).*abs(r1(range_order)),angle0,'linear',0);
        phase_a(order)=interp1(theta(range_order),phase(range_order),angle0,'linear','extrap');
        qx_a(order)=interp1(theta(range_order),qx(range_order),angle0,'linear','extrap');
        I_a2(order)=interp1(theta(range_order),I(range_order).*abs(r1(range_order)),-angle0-pi,'linear',0);
        phase_a2(order)=interp1(theta(range_order),phase(range_order),-angle0-pi,'linear','extrap');
        qx_a2(order)=-interp1(theta(range_order),qx(range_order),-angle0-pi,'linear','extrap');
    end
    Esum(num)=sum(sqrt(I_a).*exp(-1i*phase_a).*exp(-1i*2*pi*nb/lam0*qx_a.*cos(angle0)));  
    Esum2(num)=sum(sqrt(I_a2).*exp(-1i*phase_a2).*exp(-1i*2*pi*nb/lam0*qx_a2.*cos(angle0)));  % to consider the contribution of the light rays enter from the other side of the sphere
    Isum(num)=abs(Esum(num)+Esum2(num)).^2;  
    fprintf('%d/%d is completed\n',num,length(scan));
end
figure; subplot(1,2,1); plot(d,rad2deg(theta+pi/2));
hold on; subplot(1,2,2); plot(-d,-rad2deg(theta+pi/2));
figure; subplot(1,2,1); plot(d,qx);
hold on; subplot(1,2,2); plot(-d,-qx);
figure; subplot(1,2,1); plot(d,phase);
hold on; subplot(1,2,2); plot(-d,phase);
figure; subplot(1,2,1); plot(d,I.*abs(r1));
hold on; subplot(1,2,2); plot(-d,I.*abs(r1));
xx=linspace(0,30,2000);   % the x coordinate, unit: cm
zz=6.5;     % the distance between the screen and the sample, unit: cm
[X,Y]=meshgrid(xx,zz);
[thetaxy,rho]=cart2pol(X,Y);
Ir=cos(thetaxy-pi/2).^2./Y.*interp1(scan+pi,Isum,thetaxy);
x0=linspace(-10,10,1000);
y0=linspace(-10,10,1000);
[x,y]=meshgrid(x0,y0);
[phi,r]=cart2pol(x,y);
I2D=interp1(xx,Ir,r);
R=zeros(1,64);
G=zeros(1,64);
B=linspace(0,1,64);
blue_ring=[R' G' B'];
figure; imagesc(x0,y0,I2D); axis image; colormap(abs(blue_ring));