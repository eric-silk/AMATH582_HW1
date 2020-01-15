clear all; close all; clc;
load Testdata
SAMPLES = 20;
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);


if false
figure();
for j=1:SAMPLES
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    isosurface(X,Y,Z,abs(Un),0.5)
    axis([-20 20 -20 20 -20 20]), grid on
    %winput("Press anything to continue...");
end
drawnow;
end



% Plot the averaged spectrums
averaged_spectrum = zeros(n, n, n);
figure();
for j=1:SAMPLES
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftn(tmp);
    subplot(4,5,j);
    histogram(abs(fft_tmp))
    averaged_spectrum = averaged_spectrum + tmp;
end



averaged_spectrum = averaged_spectrum / SAMPLES;
figure()
histogram(abs(averaged_spectrum))
input("Enter anything to continue...")
close all
isosurface(X, Y, Z, abs(fftshift(averaged_spectrum)), 0.2)
axis([-20 20 -20 20 -20 20]); grid on, drawnow