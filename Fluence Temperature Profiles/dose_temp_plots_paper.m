clc
clear all
close all

load('Processed Profiles/high_temp_high_dose_profiles_avg.mat')

% fig = figure;
% 
% 
% % % yyaxis left 
% % plot(y,t*z(19,:),'LineWidth',1)
% % % hold on
% % % errorbar(y,z(16,:)+273,e,'.')
% % grid on
% % xlabel('Location on Sample (mm)')
% % ylabel('Fluence (m^{-2})')
% % set(gcf,'color','w');
% % set(gca,'fontsize',12);

t=1400;
y_hdht=y;
fl_hdht=t*z(19,:);
t_hdht=zt(19,:)+273.16;

clearvars -except y_hdht fl_hdht t_hdht


load('Processed Profiles/high_temp_low_dose_profiles.mat')
y_ldht=y;
fl_ldht=t*z(19,:);
t_ldht=zt(19,:)+273.16;

clearvars -except y_hdht fl_hdht t_hdht y_ldht fl_ldht t_ldht


load('Processed Profiles/low_temp_low_dose_profiles.mat')
y_ldlt=y;
fl_ldlt=t*z(19,:);
t_ldlt=zt(19,:)+273.16;

clearvars -except y_hdht fl_hdht t_hdht y_ldht fl_ldht t_ldht y_ldlt fl_ldlt t_ldlt

load('Processed Profiles/low_temp_high_dose_profiles.mat')
t=1400;
y_hdlt=y;
fl_hdlt=t*z(19,:);
t_hdlt=zt(19,:)+273.16;


close all
figure
plot(y_hdht,t_hdht,'--O',y_ldht,t_ldht,'--^',y_ldlt,t_ldlt,'--d',y_hdlt,t_hdlt,'--s','LineWidth',1)
grid on
xlabel('Radial position on Sample (mm)')
ylabel('Temperature (K)')
set(gcf,'color','w');
set(gca,'fontsize',12);
xlim([0 10])

legend('HT HD','HT LD','LT LD','LT HD')
legend('boxoff')
legend('orientation','horizontal')
legend('location','north')
% make the zeros in the fluence arrays to a vlaue and set that as ref 

k2=find(~fl_hdht);
fl_hdht(k2)=1e25;

k2=find(~fl_ldht);
fl_ldht(k2)=1e25;

k2=find(~fl_hdlt);
fl_hdlt(k2)=1e25;

k2=find(~fl_ldlt);
fl_ldlt(k2)=1e25;


figure
semilogy(y_hdht,fl_hdht,'--O',y_ldht,fl_ldht,'--^',y_ldlt,fl_ldlt,'--d',y_hdlt,fl_hdlt,'--s','LineWidth',1)
grid on
xlabel('Radial position on Sample (mm)')
ylabel('Fluence (m^{-2})')
set(gcf,'color','w');
set(gca,'fontsize',12);

legend('HT HD','HT LD','LT LD','LT HD')
legend('boxoff')
legend('orientation','horizontal')
legend('location','north')

axis([0 10 1e25 3e27])

yticks([1e25 1e26 1e27])
yticklabels({'Ref.','10^{26}','10^{27}'})
%save('temp_fluence_plots_combined_data.mat')

