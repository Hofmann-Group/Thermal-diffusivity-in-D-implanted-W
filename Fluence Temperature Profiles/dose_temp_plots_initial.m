% this code makes the linear flux profile into a 2D one 
% it also applies BC's and creates the final processed profiles of flux (z)
% and temperature (zt)
% the raw data is imported from the matlab file - temp_flux_profiles.mat
% it imports a 37x24 matrix of the temp and flux measurements for all
% samples
% select the indices such that they capture from 0 to one extreme of the
% profile 
% the description of the data for each column is given in the excel file

% for certain samples, we take the average of the profiles before and after
% irradiation, for the temperature profiles


% select the indices such that they capture from 0 to one extreme of the
% profile 

% the description of the data for each column is given in the excel file

%% need to import the profiles data 

close all
clc
clear all

load('temp_flux_profiles.mat')

% figure 
% 
% plot(profiles(1:20,1),profiles(1:20,2), 'k.','LineWidth',1)
% 
% grid on
% 
% set(gcf,'color','w');
% set(gca,'fontsize',12);
% grid on
% title('Flux Profile')

%%
% the ones with t are for the temperature profiles 

t=70*1; %time - exposed 



% the 20 shot ones, need the average of before exposure and after exposure
% for the profiles 


% LTLD
bt1=profiles(17:37,3);
bt2=profiles(17:37,4); % 

% HTHD
% bt1=profiles(16:34,23);
% bt2=profiles(16:34,24);
% 
% x=-8:0.5:8;
% y=-8:0.5:8;

% use the following 2 lines for LTHD 
% b1=(profiles(2:21,13)+profiles(1:20,9))/2;
% b2=(profiles(1:20,14)+profiles(1:20,10))/2;
% 

% this is for LTLD
b1=profiles(1:20,1);
b2=profiles(1:20,2);

% HTHD
% b1=profiles(1:20,21);
% b2=profiles(1:20,22);

x=-9:.5:9;
y=-9:.5:9;



a=zeros(length(x),length(y));
z=zeros(length(x),length(y));

% putting the data in new maticses where the values of covered region are
% adjsuted

for i=1:length(x)
    for j=1:length(y)
        
r=sqrt(x(i)^2+y(j)^2);
if r>8.5
       z(i,j)=0;
    zt(i,j)=NaN;
else
    
    k= find(abs(abs(b1)-r) < 0.9,1)
     c=b2(k);
    z(i,j)=c;
    i
    j
    
    
    kt= find(abs(abs(bt1)-r) < 0.9,1)
     ct=bt2(kt);
    zt(i,j)=ct;
    
    
end
% makes profiels go to zero for ocvered region forcefully for fluence
% for temperature it is NaN since unknown, hence we assume it is the same
% as the side since it is quite close 
if y>8.5
       z(i,j)=0;
    zt(i,j)=NaN;
end
    
    
    
    end
end

% figure
% surf(x,y,zt+273)
% title('Temperature Profile')
% xlabel('X (mm)')
% ylabel('Y (mm)')
% %zlabel('Flux (m^{-2}s^{-1})')
% zlabel('Temperature (K)')
% 
% set(gcf,'color','w');
% set(gca,'fontsize',12);
% colorbar
% 

%flux profile 
figure

surf(x,y,z)
title('Flux Profile')
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Flux (m^{-2}s^{-1})')
%zlabel('Temperature (K)')

set(gcf,'color','w');
set(gca,'fontsize',12);
colorbar

% line profile 

fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]); 

yyaxis left 
plot(y,t*z(19,:),'ko','LineWidth',1)
% hold on
% errorbar(y,z(16,:)+273,e,'.')
grid on
xlabel('Location on Sample (mm)')
ylabel('Fluence (m^{-2})')
set(gcf,'color','w');
set(gca,'fontsize',12);

d=2.0;  % distance from centre 
e=(y./y)*15;

yyaxis right 
plot(y,zt(19,:)+273,'bx','LineWidth',1)
% hold on
% errorbar(y,zt(16,:)+273,e,'.')
grid on
xlabel('Location on Sample (mm)')
ylabel('Temperature (K)')
set(gcf,'color','w');
set(gca,'fontsize',12);

legend('Fluence','Temperature')


%% saving 
% save according to the filenames that will be read by the 2D temp fleunce
% plotting code 

%save('low_temp_high_dose_profiles_2.mat','y','z','zt')


%%  additional plotting 


% %% averaging the profiles for the high dose samples 
% clc
% clear all
% close all
% load('high_temp_high_dose_profiles_before.mat','y','z','zt')
% y1=y;
% z1=z;
% zt1=zt;
% 
% load('high_temp_high_dose_profiles_after.mat','y','z','zt')
% 
% t=1400;
% 
% z_avg=0.5*(z1+z);
% zt_avg=0.5*(zt1+zt);
% clear z zt 
% z=z_avg;
% zt=zt_avg;
% 
% fig = figure;
% left_color = [0 0 0];
% right_color = [0 0 1];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]); 
% 
% yyaxis left 
% plot(y,t*z(19,:),'LineWidth',1)
% % hold on
% % errorbar(y,z(16,:)+273,e,'.')
% grid on
% xlabel('Location on Sample (mm)')
% ylabel('Fluence (m^{-2})')
% set(gcf,'color','w');
% set(gca,'fontsize',12);
% 
% d=2.0;  % distance from centre 
% e=(y./y)*15;
% 
% yyaxis right 
% plot(y,zt(19,:)+273,'LineWidth',1)
% % hold on
% % errorbar(y,zt(16,:)+273,e,'.')
% grid on
% xlabel('Location on Sample (mm)')
% ylabel('Temperature (K)')
% set(gcf,'color','w');
% set(gca,'fontsize',12);
% 
% legend('Fluence','Temperature')
% 
% save('hig_temp_high_dose_profiles_avg.mat','y','z','zt')
