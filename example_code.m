clear all; close all; clc;
load Testdata
SAMPLES = 20;
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1);
x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];
ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);


% Time Domain (I think...)
if false
for j=1:SAMPLES
    figure();
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    isosurface(X,Y,Z,abs(Un),1.25)
    axis([-20 20 -20 20 -20 20]), grid on
    %input("Press anything to continue...");
end
drawnow;
end

% Frequency Domain (again, I think...)
if false
    for j=1:SAMPLES
        figure();
        Un(:,:,:)=reshape(Undata(j,:),n,n,n);
        Un_fft=fftn(Un);
        isosurface(Kx,Ky,Kz,abs(fftshift(Un_fft)))
        %axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi])
        grid on
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
    averaged_spectrum = averaged_spectrum + fft_tmp;
end



averaged_spectrum = abs(averaged_spectrum) / SAMPLES;
figure();
title("Averaged Spectrum");
histogram(averaged_spectrum)
figure();
title("Averaged Spectrum, 100")
isosurface(Kx, Ky, Kz, fftshift(averaged_spectrum), 100)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]); 
grid on, drawnow
figure();
title("Averaged Spectrum, 125")
isosurface(Kx, Ky, Kz, fftshift(averaged_spectrum), 125)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]); 
grid on, drawnow
figure();
title("Averaged Spectrum, 150")
isosurface(Kx, Ky, Kz, fftshift(averaged_spectrum), 150)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]); 
grid on, drawnow
figure();
title("Averaged Spectrum, 175")
isosurface(Kx, Ky, Kz, fftshift(averaged_spectrum), 175)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]); 
grid on, drawnow
figure();
title("Averaged Spectrum, 200")
isosurface(Kx, Ky, Kz, fftshift(averaged_spectrum), 200)
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]);
grid on, drawnow

% Generate a 3D gaussian centered about x=1.885, y=-1.18, z=0
sigma_x = 2;
sigma_y = 2;
sigma_z = 2;
mu_x = 1.885;
mu_y = -1.18;
mu_z = 0;

coeff = 1/(2*pi*sigma_x*sigma_y*sigma_z);
exponent = (Kx-mu_x).^2/(2*sigma_x^2) + (Ky-mu_y).^2/(2*sigma_y^2) + (Kz-mu_z).^2/(2*sigma_z^2);
gauss_filter = coeff*exp(-exponent);
% normalize it
gauss_filter = gauss_filter/(max(gauss_filter(:)));
figure();
title("Filter kernel");
isosurface(Kx, Ky, Kz, gauss_filter, 0.9);
axis([-2*pi 2*pi -2*pi 2*pi -2*pi 2*pi]);
grid on, drawnow

%
figure();
for j=1:20
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftshift(fftn(tmp));
    filtered = fft_tmp.*gauss_filter;
    filtered_time = ifftn(fftshift(filtered));
    subplot(4,5,j);
    histogram(abs(filtered_time));
end
%}
figure();
title("Progress of Marble in Fluffy's innards")
xlabel("X")
ylabel("Y")
zlabel("Z")
hold on;
for j=1:20
    tmp(:,:,:) = reshape(Undata(j, :), n, n, n);
    fft_tmp = fftshift(fftn(tmp));
    filtered = fft_tmp.*gauss_filter;
    filtered_time = ifftn(fftshift(filtered));
    isosurface(X, Y, Z, abs(filtered_time), 0.3);
    %{
    hiso = patch(isosurface(X, Y, Z, abs(filtered_time), 0.3),...
                 'FaceColor', [1, 0.75, 0.65],...
                 'EdgeColor', 'none');
    isonormals(filtered_time, hiso);
    %}
    axis([-20 20 -20 20 -20 20]);
    grid on;
end