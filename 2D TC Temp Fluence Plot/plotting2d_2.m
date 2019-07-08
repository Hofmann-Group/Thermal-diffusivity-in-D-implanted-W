% this code takes in the processed profiles of tmeperature and fluence and
% the processed TGS data and created the 2D data of Thermal diffusivity
% with fluence and temperature for the high temperature samples - samples 1
% and 4
% run this code without the last two sections and it will create a file
% 'TC_fluence_temp_data_samples_2_3.mat'
% simillary plotting_2d will create the data for samples 1 and 4
% then using the last 2 sections of  plotting2d we can
% generate the main plot 



clc
clear all
close all

% q=1 is the first smaple and then q=2 it loads and does the second sample 
% first plot them and adjust the p value accordingly 
q=1;

for q=1:2
%% first loading and adjusting the TC profile
if q==1
        load('Output Data/sample3_line3_analysis_cleaned.mat','map_diffuse','std_diffuse','p')
%load('sample1_line1_analysis_cleaned.mat','map_diffuse','std_diffuse','p')
p=p-10.9;

p=p(1:end-2);
map_diffuse=map_diffuse(1:end-2); 

% for the high temp one
% p=p(1:end-5);
% map_diffuse=map_diffuse(1:end-5); 
else 
     load('Output Data/sample2_line1_analysis.mat','map_diffuse','std_diffuse','p')
    %load('sample_4_line2_analysis_cleaned.mat','map_diffuse','std_diffuse','p')
p=p-.35;
end
n=1;
n2=length(p);

%% this plot is to plot the TC profile and centre it using p above
errorbar(p(n:n2),map_diffuse(n:n2),std_diffuse(n:n2))

hold on
plot(p(n:n2),map_diffuse(n:n2), 'r','LineWidth',1)
grid on 
xlabel('Location on Sample (mm)','FontSize',14)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',12);


%% loading the flux profiles - plotting them  
if q==1
load('Processed Profiles/low_temp_low_dose_profiles.mat','y','z','zt')
t=70;
else
    
    load('Processed Profiles/low_temp_high_dose_profiles_2.mat','y','z','zt')
t=1400;
end

%% this plot is to plot the flux profile and TC 
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



%% here we take the profile data and remove the NaN in the thermal data

a=t*z(19,:); % taking the fluence data 

% taking the temperature data 
% removing the Nan value from the edges and letting the fit predict it
b=zt(19,:)+273;
[n] = find(~isnan(b));
bb=b(n);
yy=y(n);
% cftool


%% this plot is to plot the TC with the fluence profile 
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]); 
yyaxis left 
plot(p(n:n2),map_diffuse(n:n2),'LineWidth',1)
grid on 
yyaxis right 
plot(y,t*z(19,:),'LineWidth',1)

%% fitting the flux profile to gaussian 
ff = fit(y',a','gauss1');      
a1=ff.a1;
b1=ff.b1;
c1=ff.c1;

ft=fit(yy',bb','gauss1');

a2=ft.a1;
b2=ft.b1;
c2=ft.c1;


% setting beam radius 
if q==1
    
    beam_rad1=6.8;
beam_rad2=-9.1;
    
    %high temp
% beam_rad1=6.8;
% beam_rad2=-9.1;
else 
    beam_rad1=7.2;
    beam_rad2=-7.3;
end
    
for i=1:length(p)
fit1(i)=a1*exp(-((p(i)-b1)/c1)^2);
fit2(i)=a2*(exp(-((p(i)-b2)/c2)^2));

if p(i)>beam_rad1
    fit1(i)=1e25;
%     fit2=
end
if p(i)<beam_rad2
    fit1(i)=1e25;
end
end

%% plotting the fitted fleunce profiles with the measured and TC 
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]); 
yyaxis left 
plot(p(n:n2),map_diffuse(n:n2),'LineWidth',1)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',14)
yyaxis right 
plot(p,fit1,y,a)
grid on

set(gcf,'color','w');
set(gca,'fontsize',12);
xlabel('Location on Sample (mm)','FontSize',14)
ylabel('Fluence (m^{-2})','FontSize',14)
legend('TC','Fluence. Fit','Fluence. Meas.')


%% plotting the fitted temp profiles with the TC and measured temp 


fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]); 

yyaxis left 
plot(p(n:n2),map_diffuse(n:n2),'LineWidth',1)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',14)
yyaxis right 
plot(p,fit2,y,b)
grid on
set(gcf,'color','w');
set(gca,'fontsize',12);
xlabel('Location on Sample (mm)','FontSize',14)
ylabel('Temperature (K)','FontSize',14)
legend('TC','Temp. Fit','Temp. Meas.')
title('temp fit with diffuse')

%% saving the data into separate variables for the different samples 
if q==1
    map_diffuse_1=map_diffuse;
    fit1_1=fit1;
    fit2_1=fit2;
    clear map_diffuse fit1 fit2
else
    map_diffuse_2=map_diffuse;
   fit1_2=fit1;
    fit2_2=fit2;
     clear map_diffuse fit1 fit2
end

end

%% plotting the final plot 


ind=1:1:length(p);
% %% checking indices of the points that are not in the moly affected region 
% if(q==1)
% ind=find(abs(p)>4.1);
% else
%     ind=find(p<-4.1);
% end

% finding the values that are zero in the fluence and making them 1E25 for
% hte log plot 
k2=find(~fit1_1);
fit1_1(k2)=1e25;

k2=find(~fit1_2);
fit1_2(k2)=1e25;


%% plotting the low temp ones 
figure


scatter(fit1_1(ind),fit2_1(ind),120,map_diffuse_1(ind),'o','filled','LineWidth',1.5);
% set(gca,'xscale','log')
grid on
xlabel('Fluence (m^{-2})','FontSize',14)
ylabel('Temperature (K)','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',14);
colorbar
hold on
scatter(fit1_2(ind),fit2_2(ind),120,map_diffuse_2(ind),'^','filled','LineWidth',1.5);
 set(gca,'xscale','log')

hold on 

%save('TC_fluence_temp_data_samples_1_4.mat','fit1_1','fit1_2','fit2_1','fit2_2','map_diffuse_1','map_diffuse_2','p')
clear all
%% loading and plotting the high temp ones 

load('TC_fluence_temp_data_samples_1_4.mat','fit1_1','fit1_2','fit2_1','fit2_2','map_diffuse_1','map_diffuse_2','p')

% finding the values that are zero in the fluence and making them 1E25 for
% hte log plot 
k2=find(~fit1_1);
fit1_1(k2)=1e25;

k2=find(~fit1_2);
fit1_2(k2)=1e25;


%% checking indices of the points that are not in the moly affected region 
% if(q==1)
ind2=find(abs(p)>4.1);
% else
     ind3=find(p<-4.1);
% end

scatter(fit1_1(ind2),fit2_1(ind2),120,map_diffuse_1(ind2),'d','filled','LineWidth',1.5);
set(gca,'xscale','log')
grid on
xlabel('Fluence (m^{-2})','FontSize',14)
ylabel('Temperature (K)','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',14);
c=colorbar;
c.Label.String='Thermal Diffusivity (m^{2}s^{-1})';
c.FontSize=14;
c.Location='westoutside';
hold on
scatter(fit1_2(ind3),fit2_2(ind3),120,map_diffuse_2(ind3),'s','filled','LineWidth',1.5);
legend('HT HD','HT LD','LT LD','LT HD')
legend('boxoff')
legend('orientation','horizontal')
axis([1e25 .2e28 350 650])