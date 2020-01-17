% 3D 2D gaussian

x = [-128:128];
y = x;

[X, Y] = meshgrid(x,y);
sigma_x = 10;
sigma_y = 10;
mu_x = 50;
mu_y = 0;

%X, Y = meshgrid(x,y);
z = (1/2*pi*sigma_x*sigma_y)*exp(-((X-mu_x).^2/(2*sigma_x^2)+(Y-mu_y).^2/(2*sigma_y^2)));

figure()
surf(X, Y, z)