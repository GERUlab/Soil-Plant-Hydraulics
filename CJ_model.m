%Code used for calculating the hydraulic surface E(psi_leaf, psi_soil) and
%the SOL used for Carminati & Javaux, Trends in Plant Science, 2020. Code
%used for Figs.3,4,5. Plant and soil parameters to be adjusted. 
%Authors: Carminati & Javaux
%https://www.cell.com/trends/plant-science/fulltext/S1360-1385(20)30118-7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;

%------------inputs---------------
%For varying soil matric potentials hb and transpiration rates E we
%calculate the leaf water potential hleaf
hb=-[10:100:15010]';%cm bulk soil water potential in cm heads
E=3*(10^-6).*[0:1:3000]';%transpiration in cm3 s^-1

%root-plant properties
Rroot=0.15*10^7;%hPa cm^-3 s; %for the blue 0.08 % for the green 0.12
r0= 0.05; % cm root radius in cm
r2= 1; %bulk+root+rhizo soil radius in cm
L = 1200; %root length in cm - for the blue 1000 - for the green 10000

% Brooks-Corey parameters for the soil hydraulic properties
h0=-10;%this is alfa [cm^-1]
l=1.5;
ths=0.45;%not needed in the code, but useful if the decreasing water content needs to be calculated
thr=0.01;%not needed in the code, but useful if the decreasing water content needs to be calculated
k0=0.3*10^-2;%cm/s - for the blue 0.2
tau=2.5;%corresponding to 2+l*(a+2) with l exp for theta(h) and a tortusoity 

% radial rhizosphere domain
r=[r0:0.01:r2];r=r(:);
dr=(r2-r0)./(size(r)-1);
dr=dr(1);

%parameters for cavitation; if you plot the harmonic mean of 1/Rroot and k_x 
%you get a conductance that maximum is the one of root and if hx decreases the conductivity drops 
k0_x=1/Rroot; %hPa-1 cm^3 s-1;
h0_x=-20000; %in cm (or 10^-4 Mpa) approximately equal p50;
tau_x=5;
%k_x=k0_x.*(hleaf./h0_x).^-tau_x; %this would be the xylem conductance

%------------model---------------
%--> for given soil and leaf water potential we calculate the transpiration
%rate based on the conductivities. 
%here we prepare the matrix
hroot=zeros(length(hb),length(E));
hx=hroot;
hx_max=hroot;
hleaf=hroot;
Emax=hroot;
hbmat=hroot;
Emat=hroot;

%calculation of hleaf
for j1=1:length(hb)%size of psi_soil
    for j2=1:length(E)%size of E
         Emat(j1,j2)=E(j2); 
         hbmat(j1,j2)=hb(j1);
         csoil=-2*pi*r0*L* k0/(1-tau)/(h0^(-tau))/(r0/2- r0*r2^2*(log(r2)-log(r0))/(r2^2-r0^2) );
         hroot(j1,j2)=-abs(-E(j2)/csoil+hb(j1)^(1-tau)).^(1/(1-tau));            
         if (hb(j1)-hroot(j1,j2))>200000 
             break
         end
%dissipation in the root
       hx(j1,j2)=hroot(j1,j2)-Rroot.*E(j2);      % Couvreur model
%hmax - this is when the system is linear
        hx_max(j1,j2)=hbmat(j1,j2)-Rroot.*E(j2);   
%dissipation in the xylem including cavitation 
        cx=-k0_x/(1-tau_x)*(h0_x^tau_x);
        hleaf(j1,j2)=-abs(E(j2)/cx+hx(j1,j2)^(1-tau_x)).^(1/(1-tau_x));
        if (hx(j1,j2)-hleaf(j1,j2))>30000 
            break
        end           
    end
end


%Now we prepare the grid for plotting the data
[hleaf_reg,hb_reg]=meshgrid((-30000:100:0),hb); % Regular grid
E_reg=griddata(hleaf,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
Emax_reg=griddata(hx_max,hbmat,Emat,hleaf_reg,hb_reg);         % 3D Interpolates
diff_E=Emax_reg-E_reg;                                  % Diff. with Emax


%Now we calculate the SOL - where stomata close
[FX,FY]=gradient(E_reg,100,10);%FX=derivative dE/dhleaf (dhleaf is 100 in the regualr grid)
hb_reg_uni=flipud(unique(hb_reg));
hleaf_reg_uni=(unique(hleaf_reg));
gradmax=0.7;%criterion for the relative gradient - one could chose a different value (0.5, 0.8, etc... the closer it is to one, the more sensitive is the stomatal closure)
trajectory2=zeros(length(hb_reg_uni),3) ; % Trajectory between zones
   for i=1:length(hb_reg_uni)
    rel_grad(i,:)=abs(FX(i,:))./max(abs(FX(i,:)),[],2);
    pos_all=min(find(rel_grad(i,:)>=gradmax)); 
    if ~isempty(pos_all)
        pos=pos_all(end);
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(pos),E_reg(i,pos)];
    else
        trajectory2(i,:)=[hb_reg_uni(i),hleaf_reg_uni(1),E_reg(i,1)];
    end
    target=hleaf_reg(i,:);
    pos=find(target>=hb_reg_uni(i));
    E_reg(i,pos)=NaN; % impossible zone
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zon_reg=ones(size(E_reg)); % Separate zones
for i=1:length(hb_reg_uni)
    % Define zones
    target=hleaf_reg(i,:);
    pos=(target>=hb_reg_uni(i));
    zon_reg(i,pos)=0; % imposible zone
    
    target=trajectory2(i,2);
    test=hleaf_reg(i,:);
    pos=find(test<target);
    zon_reg(i,pos)=2;
end


%------------outputs/plots---------------
%Here we draw the soil isolines - E(psi_leaf) curves for constant psi_soil;
%we chose 6 lines that cover the range of soil matric potential between 0
%and -1.5 MPa
trajecthleaf1=[hleaf_reg_uni(:)*0+hb_reg_uni(1),hleaf_reg_uni(:),E_reg(1,:).'];
trajecthleaf30=[hleaf_reg_uni(:)*0+hb_reg_uni(30),hleaf_reg_uni(:),E_reg(30,:).'];
trajecthleaf50=[hleaf_reg_uni(:)*0+hb_reg_uni(48),hleaf_reg_uni(:),E_reg(48,:).'];
trajecthleaf70=[hleaf_reg_uni(:)*0+hb_reg_uni(65),hleaf_reg_uni(:),E_reg(65,:).'];
trajecthleaf100=[hleaf_reg_uni(:)*0+hb_reg_uni(100),hleaf_reg_uni(:),E_reg(100,:).'];
trajecthleaf150=[hleaf_reg_uni(:)*0+hb_reg_uni(150),hleaf_reg_uni(:),E_reg(150,:).'];

%figure from the leaf perspective
figure(1), plot(-trajecthleaf1(:,2).*10^-4,trajecthleaf1(:,3).*10^3,'k','linewidth',2)%the SOL
hold on 
plot(-trajecthleaf30(:,2).*10^-4,trajecthleaf30(:,3).*10^3,'k','linewidth',2)%isolines
plot(-trajecthleaf50(:,2).*10^-4,trajecthleaf50(:,3).*10^3,'k','linewidth',2)
plot(-trajecthleaf70(:,2).*10^-4,trajecthleaf70(:,3).*10^3,'k','linewidth',2)
plot(-trajecthleaf100(:,2).*10^-4,trajecthleaf100(:,3).*10^3,'k','linewidth',2)
plot(-trajecthleaf150(:,2).*10^-4,trajecthleaf150(:,3).*10^3,'k','linewidth',2)
plot(-trajectory2(:,2).*10^-4,trajectory2(:,3).*10^3,'r','linewidth',2)
ylabel('E [mmol m^{-2} s^{-1}]','fontsize',18,'fontname','cambria')
xlabel('-\psi_l_e_a_f [MPa]','fontsize',18,'fontname','cambria')


%Now we draw lines for isohydric curves at different psi_leaf (the value of psi_leaf is the intercept of the curves at E=0. The most negative psi_leaf here is -3MPa)
trajecthsoil200=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(281),E_reg(:,281)];%psi_leaf=-2.0MPa
trajecthsoil500=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(251),E_reg(:,251)];%psi_leaf=-1.5MPa
trajecthsoil750=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(226),E_reg(:,226)];%psi_leaf=-4.0MPa
trajecthsoil1000=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(201),E_reg(:,201)];%psi_leaf=-1.0MPa
trajecthsoil1500=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(151),E_reg(:,151)];%psi_leaf=-0.5MPa
trajecthsoil3000=[hb_reg_uni(:),hb_reg_uni(:)*0+hleaf_reg_uni(1),E_reg(:,1)];%psi_leaf=-0.3MPa


%fig from soil perspective 
figure(2),%,plot(Pr(:,1).*10^-4,Tr(:,1),'o',Pr(:,2).*10^-4,Tr(:,2),'o',Pr(:,3).*10^-4,Tr(:,3),'o',Pr(:,4).*10^-4,Tr(:,4),'o',Pr(:,5).*10^-4,Tr(:,5),'o',Pr(:,6).*10^-4,Tr(:,6),'o',Pr(:,7).*10^-4,Tr(:,7),'o',Pr(:,8).*10^-4,Tr(:,8),'o',Pr(:,9).*10^-4,Tr(:,9),'o')
hold on
plot(-trajectory2(:,1).*10^-4,trajectory2(:,3),'r','linewidth',2)%The SOL
plot(-trajecthsoil750(:,1).*10^-4,trajecthsoil750(:,3),'k','linewidth',2)%isolines
plot(-trajecthsoil200(:,1).*10^-4,trajecthsoil200(:,3),'k','linewidth',2)
plot(-trajecthsoil500(:,1).*10^-4,trajecthsoil500(:,3),'k','linewidth',2)
plot(-trajecthsoil1000(:,1).*10^-4,trajecthsoil1000(:,3),'k','linewidth',2)
plot(-trajecthsoil1500(:,1).*10^-4,trajecthsoil1500(:,3),'k','linewidth',2)
plot(-trajecthsoil3000(:,1).*10^-4,trajecthsoil3000(:,3),'k','linewidth',2)
ylabel('E [mmol m^{-2} s^{-1}]','fontsize',18,'fontname','cambria')
xlabel('-\psi_s_o_i_l [MPa]','fontsize',18,'fontname','cambria')

%calculation of psi_s50, psi_gs50, psi_x50 and psi_H50 - note that psi_gs50
%and psi_s50 are function of the maximum transpiration rate max(E). These
%definitions were used to prodice Fig.5
i=find(trajectory2(:,3)<= (trajectory2(1,3)./2+trajectory2(1,3)./20) & trajectory2(:,3)>= (trajectory2(1,3)./2-trajectory2(1,3)./20));
i=round(median(i));
h_gs50=trajectory2(i,2)
h_s50=-abs(0.5*max(E)./csoil).^(1/(1-tau))
h_x50=h0_x
h_H50=max(h_s50,h_x50)


%3D surface with SOL
figure
mapsurf=[0 0 0  ; 35 168 48 ; 127 91 14 ]./256;
surface(hb_reg.*10^-4,hleaf_reg.*10^-4,E_reg,zon_reg);hold on
colormap(gca,mapsurf);
view([-120,25])
shading flat;grid on
plot3(trajectory2(:,1).*10^-4,trajectory2(:,2).*10^-4,trajectory2(:,3),'r','linewidth',2)
hold off
zlabel('E [mmol m^{-2} s^{-1}]','fontsize',18,'fontname','cambria')
xlabel('-\psi_s_o_i_l [MPa]','fontsize',18,'fontname','cambria')
ylabel('-\psi_l_e_a_f [MPa]','fontsize',18,'fontname','cambria')
