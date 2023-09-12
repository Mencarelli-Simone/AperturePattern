clear;
clc;
close all
fclose all;

L = 1.5;    % length of antenna
W = 0.3;    % width of antenna

lambda = 0.03; % 10 GHz
k = 2*pi/lambda;

dx = lambda / 5;       % size of the integration bin
dy = lambda / 5;

Nx = W / dx;
Ny = L / dy;

antenna_distorted = zeros(Nx,Ny,3);

% assume a nominal curved antenna shape
for ii = 1:Nx
    for jj = 1:Ny
        antenna_distorted(ii,jj,1) = -W/2 + (ii-1)*dx;
        antenna_distorted(ii,jj,2) = -L/2 + (jj-1)*dy;
        antenna_distorted(ii,jj,3) = 0; %-0.5*lambda * (cos(jj/Ny * pi/2)-1);
    end
end

plot3(antenna_distorted(:,:,1),antenna_distorted(:,:,2),antenna_distorted(:,:,3),'b.')
axis image
drawnow

phi = [0 pi/2]; 
theta = linspace(deg2rad(-90),deg2rad(90),200);
R = 1e5;        % 100 km

f = zeros(length(theta),length(phi));

for tt = 1:length(theta)
    tt
    for pp = 1:length(phi)
       pp
        [x_ff,y_ff,z_ff] = sph2cart(phi(pp),(pi/2 - theta(tt)),R);

        for ii = 1:Nx 
            for jj = 1:Ny
                r = norm(squeeze(antenna_distorted(ii,jj,:)) - [x_ff;y_ff;z_ff]);

                f(tt,pp) = f(tt,pp) + exp(-1i*k*r);
            end
        end

    end
end

%%
figure;
plot(rad2deg(theta),20*log10(abs(f(:,1)/(Nx*Ny))));
hold on;
plot(rad2deg(theta),20*log10(abs(sinc(W*sin(theta)/lambda))),'--');     % theoretical
ylim([-60 0])
xlabel('Elevation')
ylabel('Normalised pattern (Azimuth = 90)')
legend('Aperture integration','Theoretical')

figure

plot(rad2deg(theta),20*log10(abs(f(:,2)/(Nx*Ny))));
hold on;
plot(rad2deg(theta),20*log10(abs(sinc(L*sin(theta)/lambda))),'--');     % theoretical
ylim([-60 0])
xlabel('Elevation')
ylabel('Normalised pattern (Azimuth = 0)')
legend('Aperture integration','Theoretical')
