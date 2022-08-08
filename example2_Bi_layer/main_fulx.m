clc; clear; 

x_mesh = 300;
x_min = 0;
x_max = 2;

y_min = 0;
y_max = 0.5;

max_flux = 0.15;
min_flux = 0.03; 

spike_loc = (x_min + x_max)/4;
sharpness = 50;

%%
x = linspace(x_min,x_max,x_mesh);
y = - sharpness * ( x - spike_loc ).^2 + max_flux;
rang = sqrt((max_flux - min_flux) / sharpness);
y(x < (spike_loc - rang)) = min_flux;
y(x > (spike_loc + rang)) = min_flux;

figure(1); clf
plot(x,y,'-')
axis equal
axis([min(x),max(x),y_min,y_max])
grid on
