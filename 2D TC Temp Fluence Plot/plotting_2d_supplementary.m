clc
clear all
close all

% this code will generate figure 4 as well as the supplementary figure S10,
% with error bars. It requires several data files to run, that are provided
% in the home folder, as well as in the Processed Profiles and Output Data
% folders 
% q=1 is the first smaple and then q=2 it loads and does the second sample 
% first plot them and adjust the p value accordingly 
q=1;

for q=1:2
%% first loading and adjusting the TC profile
if q==1
load('Output Data/Data/sample1_line1_analysis_cleaned.mat','map_diffuse','std_diffuse','p')
% p=p-11.5;
p=p-10.25;   % after review temp issue


p=p(1:end-5);
map_diffuse=map_diffuse(1:end-5);
else 
    load('Output Data/Data/sample_4_line2_analysis_cleaned.mat','map_diffuse','std_diffuse','p')
p=p-10;
end
n=1;
n2=length(p);

%% this plot is to plot the TC profile and centre it using p above
figure
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
load('Processed Profiles/high_temp_high_dose_profiles_avg.mat','y','z','zt')
t=1400;
else
    
    load('Processed Profiles/high_temp_low_dose_profiles.mat','y','z','zt')
t=70;
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



%% here we take the profile data and remove the NaN in the thermal data,
%remove Moly regions in thermal data and take fleunce data

a=t*z(19,:); % taking the fluence data 

% taking the temperature data 
% removing the Nan value from the edges and letting the fit predict it -
% also removing moly affected regions - new 
% 
% b=zt(19,:)+273;
% [n] = find(~isnan(b));
% bb=b(n);
% yy=y(n);

%new 
b=zt(19,:)+273;

ind_m=find(abs(y)>5.1);
bb2=b(ind_m);
yy2=y(ind_m);

[n] = find(~isnan(bb2));
bb=bb2(n);
yy=yy2(n);
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

%% fitting the flux and temp profiles to gaussian 
ff = fit(y',a','gauss1');      
a1=ff.a1;
b1=ff.b1;
c1=ff.c1;

% fitting temp profile to gaussian
ft=fit(yy',bb','gauss1');

a2=ft.a1;
b2=ft.b1;
c2=ft.c1;


% setting beam radius 
if q==1
% beam_rad1=6.8;
% beam_rad2=-9.1;

beam_rad1=7.9;
beam_rad2=-7.9;

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
    fit1_1=fit1;   % fluence 
    fit2_1=fit2;   % temperature 
    clear map_diffuse fit1 fit2
else
    map_diffuse_2=map_diffuse;
   fit1_2=fit1;    % fluence 
    fit2_2=fit2;   % temperature 
     clear map_diffuse fit1 fit2
end

end

%% plotting the final plot 

%% checking indices of the points that are not in the moly affected region 
if(q==1)
ind=find(abs(p)>4.1);
else
    ind=find(p<-4.1);
end

%% plotting the high temp ones 
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
xticks([1e25 1e26 1e27])
xticklabels({'Ref.','10^{26}','10^{27}'})
hold on 


%saving the variables for the suppl plot 
hind=ind;
hfit1_1=fit1_1;

hfit2_1=fit2_1;

hmap_diffuse_1=map_diffuse_1;

hfit1_2=fit1_2;

hfit2_2=fit2_2;

hmap_diffuse_2=map_diffuse_2;


%save('TC_fluence_temp_data_samples_1_4.mat','fit1_1','fit1_2','fit2_1','fit2_2','map_diffuse_1','map_diffuse_2','p')
clear fit1_1 fit2_1 map_diffuse_1 fit1_2 fit2_2 map_diffuse_2
%% loading and plotting the low temp ones 

load('TC_fluence_temp_data_samples_2_3.mat','fit1_1','fit1_2','fit2_1','fit2_2','map_diffuse_1','map_diffuse_2','p')

% finding the values that are zero in the fluence and making them 1E25 for
% hte log plot 
k2=find(~fit1_1);
fit1_1(k2)=1e25;

k2=find(~fit1_2);
fit1_2(k2)=1e25;


ind=1:1:length(p);
% ind1=1:1:length(p)-40;
% % ind1=1:1:length(p);
% % 

% removes the central drop points in LTHD
ind1=[1:1:40,46:1:length(p)];


scatter(fit1_1(ind),fit2_1(ind),120,map_diffuse_1(ind),'d','filled','LineWidth',1.5);
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
scatter(fit1_2(ind1),fit2_2(ind1),120,map_diffuse_2(ind1),'s','filled','LineWidth',1.5);
% legend('HT HD','HT LD','LT LD','LT HD')
% legend('boxoff')
% legend('orientation','horizontal')
axis([1e25 .2e28 350 650])


%% plot for supplementary - new figure 


% getting the error values 

    load('Output Data/Data/sample3_line3_analysis_cleaned.mat','std_diffuse')
LTLD_std_diffuse=std_diffuse;
clear std_diffuse

    load('Output Data/Data/sample2_line1_analysis.mat','std_diffuse')
LTHD_std_diffuse=std_diffuse;
clear std_diffuse

    load('Output Data/Data/sample1_line1_analysis_cleaned.mat','std_diffuse')
HTHD_std_diffuse=std_diffuse;
clear std_diffuse

    load('Output Data/Data/sample_4_line2_analysis_cleaned.mat','std_diffuse')
HTLD_std_diffuse=std_diffuse;
clear std_diffuse




% 2d figure starts here 

figure 



%low temp ones first 

load('TC_fluence_temp_data_samples_2_3.mat','fit1_1','fit1_2','fit2_1','fit2_2','map_diffuse_1','map_diffuse_2','p')




% getting the errorbars
h1=errorbar(fit1_1(ind),map_diffuse_1(ind),LTLD_std_diffuse(ind),'bd','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on





% finding the values that are zero in the fluence and making them 1E25 for
% hte log plot 
k2=find(~fit1_1);
fit1_1(k2)=1e25;

k2=find(~fit1_2);
fit1_2(k2)=1e25;


ind=1:1:length(p);
% ind1=1:1:length(p)-40;
% % ind1=1:1:length(p);
% % 

% removes the central drop points in LTHD
ind1=[1:1:42,48:1:length(p)];


scatter(fit1_1(ind),map_diffuse_1(ind),60,fit2_1(ind),'d','filled','LineWidth',1.5);
set(gca,'xscale','log')
grid on
xlabel('Fluence (m^{-2})','FontSize',14)
 ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',14);
c=colorbar;
c.Label.String='Temperature (K)';
c.FontSize=14;
c.Location='eastoutside';

% making it not come in the legend 
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


hold on


% getting the errorbars
h2=errorbar(fit1_2(ind1),map_diffuse_2(ind1),LTHD_std_diffuse(ind1),'bs','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on

% making it not come in the legend 
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



scatter(fit1_2(ind1),map_diffuse_2(ind1),60,fit2_2(ind1),'s','filled','LineWidth',1.5);
% legend('HT HD','HT LD','LT LD','LT HD')
% legend('boxoff')
% legend('orientation','horizontal')
axis([1e25 .2e28 4.8e-5 7e-5])

hold on


% HTLD
% getting the errorbars
h4=errorbar(hfit1_2(hind),hmap_diffuse_2(hind),HTLD_std_diffuse(hind),'r^','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on

% making it not come in the legend 
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


scatter(hfit1_2(hind),hmap_diffuse_2(hind),60,hfit2_2(hind),'^','filled','LineWidth',1.5);
 set(gca,'xscale','log')
xticks([1e25 1e26 1e27])
xticklabels({'Ref.','10^{26}','10^{27}'})
hold on 


% getting the errorbars
h3=errorbar(hfit1_1(hind),hmap_diffuse_1(hind),HTHD_std_diffuse(hind),'ro','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on

% making it not come in the legend 
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



scatter(hfit1_1(hind),hmap_diffuse_1(hind),60,hfit2_1(hind),'o','filled','LineWidth',1.5);
 set(gca,'xscale','log')
grid on
xlabel('Fluence (m^{-2})','FontSize',14)
%  ylabel('Temperature (K)','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',14);




legend('LTLD','LTHD','HTLD','HTHD','Location','northeast','Orientation','horizontal')

hold off 

