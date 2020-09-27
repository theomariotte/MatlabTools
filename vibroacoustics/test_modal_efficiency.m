% test modal efficiency
clear; clc; close all;

%% gloabl variables

global Lx;
global Ly;
global rho_0;
global c_sound;
global dx
global dy;

%% Parameters

% plate dimensions
Lx = 0.414;
Ly = 0.314;



% Spatial resolution
dx = 0.005;
dy = 0.005;

% speed of sound [m/s]
c_sound = 342;
% density (20°C) [kg/m^3]
rho_0 = 1.292;

%% TEST 01 : COMPUTE MODAL EFFICIENCY AS A FUNCTION OF MODE ORDERS

% Maximum orders...
% Over X
n_max = 10;
% Over Y
m_max = 10;

% acoustical frequency [Hz]
f = 50;
w = 2*pi*f;
sigma_mn = zeros(n_max,m_max);

for n = 1 : n_max
    for m = 1 : m_max



        [sigma_mn(n,m),kx,ky] = getModalEfficiency(n,m,w);

        %%%% Figures %%%%
        % eigen shape
%         figure('Name','Eigen shape');
%         im = imagesc(x,y,phi_mn);
%         xlabel('x')
%         ylabel('y')
%         ax = gca;
%         set(ax,'ydir','normal');

%         % Overall reponse
%         figure('Name','TF 2D All');
%         imagesc(kx,ky,PHI_mn2);
%         xlabel('k_x');
%         ylabel('k_y');
%         colorbar;

%         % "real" spectrum
%         PHI_mn2_real = PHI_mn2(1:floor(Nx/2)+1, 1:floor(Ny/2)+1);
%         figure('Name','TF 2D Real');
%         im = imagesc(kx(1:floor(Nx/2)+1),ky(1:floor(Ny/2)+1),PHI_mn2_real);
%         xlabel('k_x');
%         ylabel('k_y');
%         shading('interp')
%         ax = gca;
%         ax.CLim = [0 100];
%         set(ax,'ydir','normal');
%         colorbar;
    end
end

n_vec = 1 : n_max;
m_vec = 1 : m_max;

figure
pcolor(n_vec,m_vec,sigma_mn);
xlabel('n')
ylabel('m')

figure
plot(n_vec,sigma_mn(:,1),'kx','linewidth',3)
xlabel('n')
ylabel('\sigma_{mn}')
set(gca,'fontsize',14);
grid on

figure
plot(m_vec,sigma_mn(1,:),'kx','linewidth',3)
xlabel('n')
ylabel('\sigma_{mn}')
set(gca,'fontsize',14);
grid on


%% TEST 02 : COMPUTE MODAL EFFICIENCY AS A FUNCTION OF ACOUSTICAL FREQUENCY

fmin = 300;
fmax = 4e3;
df = 1;
w_vec = 2*pi*( fmin : df : fmax )';
num_freq = length(w_vec);

% order over x
n = 1;
% order over y
m = 1;

sigma_mn2 = zeros(num_freq,1);
for ii = 1 : num_freq

    [sigma_mn2(ii),kx,ky] = getModalEfficiency(n,m,w_vec(ii));

end

figure
plot(w_vec/2/pi,sigma_mn2,'k','linewidth',3)
xlabel('Frequency [Hz]')
ylabel('\sigma_{mn}')
set(gca,'fontsize',14);
grid on



%% External functions

function [sigma_mn,kx,ky] = getModalEfficiency(n,m,w)
        
    global rho_0;
    global c_sound; 
    global dx
    global dy;
    
    [phi_mn,x,y] = getEigenShape(n,m);
    
    Nx = length(x);
    Ny = length(y);

    kx = linspace(0,1/dx,Nx);
    ky = linspace(0,1/dy,Ny);
    
    PHI_mn = fft2(phi_mn);
    PHI_mn2 = abs(PHI_mn).^2;


    [KX,KY] = meshgrid(kx,ky);
    k = w/c_sound;
    k2 = k*k;
    M_tmp = PHI_mn2 ./ sqrt(k2 - KX.^2 - KY.^2);

    sigma_mn = 1/2/pi * real(w * rho_0 * sum(sum(M_tmp)));
   
end

function [phi_mn,x,y] = getEigenShape(n,m)
    global Lx;
    global Ly;
    global dx
    global dy;

    x = (0 : dx : Lx)';
    y = (0 : dy : Ly)';
    
    [X,Y] = meshgrid(x,y);
    
    invLx = 1/Lx;
    invLy = 1/Ly;
    
    phi_mn = sin(n*pi*invLx*X) .* sin(m*pi*invLy * Y);
end

