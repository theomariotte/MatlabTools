%%                      Version 11  - 14/10/2020
%             Ecrit par Téo Dubois - Apprenti ingénieur NVH, LV
%                            Octobre 2020
%--------------------------------------------------------------------------
% Ce code permet de modéliser une plaque simple et de regarder son
% comportement modale et de voir son rayonnement. Dans une première approche
% nous regarderons le cas d'une plaque simple rectangulaire, simplement
% appuyé sur son pourtour. Les paramètres seront choisis par l'utilisateur.
% 

% Travaux en cours sur cette version :
% - Tentative de mise en place du calcul de l'efficacité modale ok ?
% - Résoudre le pb avec les nombre d'onde -->  C.Pezerat ?

% Sur interface graphique : Corriger bugs qui allonge le temps de calcul
%--------------------------------------------------------------------------
%%                            Nettoyage
%--------------------------------------------------------------------------
close all;
clear variables; 
clc;
%%                           Setup MATLAB
%--------------------------------------------------------------------------
matlab_home = 'C:\Users\tdubois1\Documents\Teo_Travail_Local\01_Projet\-11-Modele_Plaque\MATLAB\';
fig_path = [matlab_home 'Figures\'];
data_path = [matlab_home 'Data\'];
%%                 Paramètres intrasèque de la plaque
%--------------------------------------------------------------------------
% --- Dimensions ---
Lx=140e-3;    % longueur [m]
Ly=70e-3; % Largeur [m]
Lz = 2e-3;     % Epaisseur [m] vérifier ok

% --- Paramètres intrasèque de la plaque ---
rho = 1730;   % Masse volumique [kg/m3]
E = 18e9;      % Module Young
nu = 0.37;    % Coefficient de poisson

% --- Fréquence d'excitation ---
f = 1700;     
f_p=10:0.1:12000; % Vecteur fréquence

% --- Paramètres du milieu fluide ---
rho_air = 1.22; % Masse volumique du fluide, ici l'air [kg/m^3]
c = 340; % Vitesse du son dans l'air [m/s²]

% --- Force injectée ---
Force = 1; % [N]
xf0 = Lx/2; % [m]
yf0 = Ly/2; % [m]

% --- Ordre des modes de vibration
m_v = 1;
m_vv=1:1:3;
n_vv=m_vv;
n_v = 1;

% Theta
% theta =linspace(0,pi*2,100000);
theta=1:0.1:2*pi;

% Phi
% phi = linspace(0,2*pi,100000);
phi=1:0.1:90;

%%              Définition de la géométrie de la plaque
%--------------------------------------------------------------------------
ds=1e-3; % Pas de discrétisation

x=0:ds:Lx; % Vecteur x de la plaque discrétisée
Nx=length(x);

y=0:ds:Ly;
Ny=length(y);

z=0:ds:Lz;
Nz=length(z);

%%                         Calculs divers
%--------------------------------------------------------------------------
% Rigidité de flexion
D = E*Lz^3/(12*(1-nu^2));

% Pulsation du Déplacement transversal
omega = f*2*pi; % Omega excitation

% Fréquence critique de la plaque
omega_critique = c^2*sqrt((rho*Lz)/D);
% lambda_critique=f/omega_critique;
freq_critique=(c^2/(2*pi))*sqrt(rho*Lz/D);
lambda_critique=c/freq_critique;
% Relation de dispertion
kf = omega*sqrt(rho*Lz/D); % Il s'agit de kf^2 = kx^2+ky^2

% Nombre d'ondee
k = omega/c;

% Longueur d'onde
lambda = c/f;

% Fréquence naturelle du mode m
%Page 44 du Jc Pascal MV2 VIBROACOUSTIQUE DES STRUCURE PLANES

% Construction de kx et ky
% kx = 0:2:100;
% ky=kx;
kx= 0:0.1*pi/Lx:100;
% ky= 0:0.1*pi/Ly:100;
ky=kx;

% Masse surfacique
rho_s=rho*Lz; 

%%              Calucls des déformées propres de la plaque 
%--------------------------------------------------------------------------
phi_mn = zeros(Nx,Ny);
for ii = 1 : Nx
    for jj = 1 : Ny
        phi_mn(ii,jj) = sin((m_v*pi*x(ii))/Lx) * ...
            sin((n_v*pi*y(jj))/Ly);
    end
end
clear ii;clear jj;

%%            Calcul des fréquences propres des modes m,n
%--------------------------------------------------------------------------
omega_pq=sqrt(D/(rho*Lz))*((m_v*pi/Lx)^2+(n_v*pi/Ly)^2);
for ii=1:length(m_vv)
    for jj=1:length(n_vv)
        omega_pq2(ii,jj)=sqrt(D/(rho*Lz))*((m_vv(ii)*pi/Lx)^2+(n_vv(jj)*pi/Ly)^2);
    end 
end
frequence_pq=omega_pq2./(2*pi);
%%      Calcul du spectre du nombre d'onde de la déformé du mode M
%--------------------------------------------------------------------------
% FONCTIONNE
% for ii=1:length(kx)
%     for jj=1:length(ky)
%         fun1= @(x_) sin(m_v*pi*x_).*exp(1i*(kx(ii)*x_));
%         term1(ii,jj)=integral(fun1,-Lx,Lx);
%         fun2= @(y_) sin(n_v*pi*y_).*exp(1i*(ky(jj)*y_));
%         term2(ii,jj)=integral(fun2,-Ly,Ly);     
%     end
% end
% clear ii; clear jj;
% phi_chapeau=term1*term2;
% phi_chapeau2=abs(phi_chapeau).^2;

for pp=1:length(m_vv)
    for qq=1:length(n_vv)
        for ii=1:length(kx)
            for jj=1:length(ky)
                fun1= @(x_) sin(m_vv(pp)*pi*x_).*exp(1i*(kx(ii)*x_));
                term1(ii,jj)=integral(fun1,-Lx,Lx);
                fun2= @(y_) sin(n_vv(qq)*pi*y_).*exp(1i*(ky(jj)*y_));
                term2(ii,jj)=integral(fun2,-Ly,Ly);
            end
        end
        phi_chapeau{pp,qq}=term1.*term2;
        phi_chapeau2{pp,qq}=abs(phi_chapeau{pp,qq}).^2;
    end
end
clear ii; clear jj;clear pp;clear qq;



%%        Recherche du maximum pour faire varier les modes
%--------------------------------------------------------------------------
for ii=1:length(m_vv)
    for jj=1:length(n_vv)
        maximum=max(max(phi_chapeau2{ii,jj}));
        [te1,te2]=find(phi_chapeau2{ii,jj}==maximum);
        toplot{ii,jj}=[te1,te2];
    end
end

%%                    Calcul de la norme
%--------------------------------------------------------------------------
Norme=Lx*Ly/4;
%%                    Calcul force généralisée
% --------------------------------------------------------------------------
Force_g=zeros(length(m_v),length(n_v));
for i=1:length(m_vv)
    for j=1:length(n_vv)
        Force_g=Force*sin(m_vv(i)*pi*xf0/Lx)*sin(n_vv(j)*pi*yf0/Ly);
    end
end
clear i; clear j;

%%               Calcul du déplacement de la plaque
%--------------------------------------------------------------------------
temp=zeros(Nx,Ny);
w=zeros(Nx,Ny);
w_tot={temp};
% for ii=1:length(m_v)
%     for jj=1:length(n_v)
        for i=1:Nx
            for j=1:Ny
                var_1=phi_mn;
                w(i,j)=(Force_g*...
                    var_1(i,j)/...
                    (Norme*rho*Lz*(omega_pq^2-omega^2)));
                temp(i,j)=w(i,j);
                w(i,j)=w(i,j)+temp(i,j); 
            end
        end
        w_tot=temp;
        w_tot=w_tot.*(4*omega/(rho_s*Lx*Ly));
%     end
% end
clear i; clear j; clear ii; clear jj;

%%                  Calcul de l'efficacité modale
%            !!!!!!!!!!!!!!!! EN TEST !!!!!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
clc;
% for ii=1:length(kx)
%     for jj=1:length(ky)
%         fun3= @(x_) real(phi_chapeau2{1,1}(1,1))./...
%                        ((k^2-kx(ii)^2-ky(jj)^2)*2*pi*2*pi);
%         term3(ii,jj)=integral(fun3,0,Lx);
%     end
% end
% % sigma=real(omega*rho*term3);
% clear ii; clear jj;
for pp=1:length(m_vv)
    for qq=1:length(n_vv)
        for ii=1:length(kx)
            for jj=1:length(ky)
                fun3=@(x_,y_) abs(sin(m_vv(pp)*pi*x_).*exp(1i*(kx(ii)*x_))...
                    .*sin(n_vv(qq)*pi*y_).*exp(1i*(ky(jj)*y_))).^2./...
                    (sqrt(k^2-kx(ii)^2-ky(jj)^2)*2*pi*2*pi);
                term3(ii,jj)=integral2(fun3,0,Lx,0,Ly);
                

            end
        end
%         phi_chapeau3{pp,qq}=term1*term2;
%         phi_chapeau3{pp,qq}=real(term3)*omega*rho_air;
        phi_chapeau3{pp,qq}=term3.*omega.*rho_air;
    end
end
clear i; clear j; clear ii; clear jj;
%% Calcul de la moyenne pour avoir le coef ??
%            !!!!!!!!!!!!!!!! EN TEST !!!!!!!!!!!!!!!!!!!

for ii=1:length(phi_chapeau3)
    for jj=1:length(phi_chapeau3)
        moy_phi3(ii,jj)=mean(abs(phi_chapeau3{ii,jj}),'All')^2;
    end
end

%% Facteur de rayonnement approché pour une plaque rectangulaire
%--------------------------------------------------------------------------
% g1= @(alpha_) (4/pi^4)*((1-2*alpha_^2)/(alpha_*sqrt(1-alpha_^2)));
% g2= @(alpha_) (1/(4*pi^2))*(((1-alpha_^2)*log((1+alpha_)/(1-alpha_))+2*alpha_)/...
% (1-alpha_^2)^(3/2));
r=Lx/Ly;    
for ii=1:length (f_p)  
    alpha(ii)=sqrt(f_p(ii)/freq_critique);
%     f_p(ii)

    g1=  (4/pi^4)*((1-2*alpha(ii)^2)/(alpha(ii)*sqrt(1-alpha(ii)^2)));
    g2=  (1/(4*pi^2))*(((1-alpha(ii)^2)*log((1+alpha(ii))/(1-alpha(ii)))+2*alpha(ii))/...
        (1-alpha(ii)^2)^(3/2));
    if f_p(ii)<floor(freq_critique)
        if f_p(ii)<0.5*freq_critique
            sigma(ii)=2*r*(lambda_critique/Lx)^2*g1+2*(1+r)*(lambda_critique/Lx)...
            *g2;
        else
            sigma(ii)=2*(1+r)*(lambda_critique/Lx)*g2;
        end
%     elseif f_p(ii)>=floor(freq_critique-20) && f_p(ii)<=floor(freq_critique+20)
%         sigma(ii) = sqrt(Lx/lambda_critique)*(1+1/sqrt(r));
    else
        sigma(ii)=(1-freq_critique/f_p(ii))^(-1/2);
    end
end
% sigma=sigma./max(sigma);
semilogx(f_p,10*log10(sigma),'Color','k','LineWidth',2)
grid on;
title('Facteur de rayonnement approché','FontSize',14);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Modal efficiency [dB]','FontSize',14);
    
%%                               Affichage 
%--------------------------------------------------------------------------
% Permet d'afficher uniquement les modes de plaque en 2D. En effet ici, il
% n'y a pas d'amplitude modale

%% Fenêtre 1
h1 = figure('Name','Déformés modales');
surf(y,x,phi_mn)
xlabel('X [m]')
ylabel('Y [m]')
shading('interp')
title('Mode de plaque');

%% Fenetre 2
h2 = figure('Name','');
plot(kx,phi_chapeau2{1,1}(1,:));
grid on;

%% Fenetre 3
h3 = figure('Name','');
surf(abs(phi_chapeau2{1,1}))
shading('interp')

%% Fenetre 4
h4 = figure('Name',"Vibroacoustique d'une plaque simple");
pcolor(abs(phi_chapeau2{1,1}))
shading('interp')
xlabel('K_x','FontSize',12);
ylabel('K_y','FontSize',12);
title(['Spectre de nombre d"onde pour une plaque pour les modes M = ',...
    num2str(m_v), ' & N = ', num2str(n_v) ' | Fréquence : ',num2str(f),...
    ' Hz'],'FontSize',12);
hold on;
grid on; 

% for ii=1:length(toplot)
%     for jj=1:length(toplot)
%         %         x_toplot=toplot{ii,jj}(1,1)
%         %         y_toplot=toplot{ii,jj}(1,2)
%         %         phi_chapeau2{ii,jj}(x_toplot,y_toplot)
%         %         plot(phi_chapeau2{ii,jj}(x_toplot,y_toplot),'r*')
%         x_toplot=toplot{ii,jj}(1,1);
%         y_toplot=toplot{ii,jj}(1,2);
%         plot(y_toplot,x_toplot,'r*')
%     end
% end

aff1=plot(k*cos(phi),k*sin(phi),'LineWidth',2,'Color','k');
aff2=plot(sqrt(kf)*cos(phi),sqrt(kf)*sin(phi),'LineWidth',2,'Color','r');
axis('equal')
xlim([1 length(kx)])
ylim([1 length(ky)])
lgd=legend([aff1 aff2], ["Nombre d'onde acoustique",...
    "Nombre d'onde vibratoire"],'FontSize',12);
title(lgd,"Légende")
x0_pos_fig=100;
y0_pos_fig=50;
width=720 ;
height=720;
set(gcf,'position',[x0_pos_fig,y0_pos_fig,width,height])
% text(10,90,['Fréquence : ' num2str(f),'LineWidth',12]);
% print(h4,[fig_path sprintf('Mode_plaque_%d_%d',m_v,n_v)],'-dpng','-r400')

%% Figure 5
h5 = figure('Name',"Vibroacoustique d'une plaque simple");
% pcolor(abs(phi_chapeau2{10,10}))
% shading('interp')
xlabel('K_x','FontSize',12);
ylabel('K_y','FontSize',12);
title(['Spectre de nombre d"onde pour une plaque | Fréquence : ',...
    num2str(f),' Hz'],'FontSize',12);
hold on;
grid on; 
aff1=plot(k*cos(phi),k*sin(phi),'LineWidth',2,'Color','k');
aff2=plot(sqrt(kf)*cos(phi),sqrt(kf)*sin(phi),'LineWidth',2,'Color','r');
for ii=1:length(toplot)
    for jj=1:length(toplot)
        %         x_toplot=toplot{ii,jj}(1,1)
        %         y_toplot=toplot{ii,jj}(1,2)
        %         phi_chapeau2{ii,jj}(x_toplot,y_toplot)
        %         plot(phi_chapeau2{ii,jj}(x_toplot,y_toplot),'r*')
        x_toplot=toplot{ii,jj}(1,1);
        y_toplot=toplot{ii,jj}(1,2);
        plot(y_toplot,x_toplot,'Marker','o','Color','black','MarkerSize',15)
    end
end
axis('equal')
xlim([1 length(kx)])
ylim([1 length(ky)])
lgd=legend([aff1 aff2], ["Nombre d'onde acoustique",...
    "Nombre d'onde vibratoire"],'FontSize',12);
title(lgd,"Légende")
x0_pos_fig=100;
y0_pos_fig=50;
width=720 ;
height=720;
set(gcf,'position',[x0_pos_fig,y0_pos_fig,width,height])







