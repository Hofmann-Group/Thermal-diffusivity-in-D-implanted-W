%% Import fixed TG data files
 clc
 clear all;
  close all;


%loading TGS data

load('TGS data/sample_3_line2.mat','yy1','yy','xx1','xx','p')

% main code uses xx yy as 2 d matrices , nmow we have 3d since locaiton is
% also changing 

%% data operations

%loading the 3D data into new matrices
% xxx1=xx1;
% xxx=xx;
% yyy=yy;
% yyy1=yy1;
% p=p;

% for sample 1 line 1 
xxx1=xx1(:,:,1:end);
xxx=xx(:,:,1:end);
yyy=yy(:,:,1:end);
yyy1=yy1(:,:,1:end);
p=p;


% plots 
splot=0;  % plots the initial plots 
iiiplot=0;  % pltos for each startin point with the main fit 
nnnplot=0 ;  % plots for each trace - with the master fit 

saving=0;

clear xx xx1 yy1 yy

% variable for saving the diffusivity with location
map_diffuse=zeros(1,length(p));

% extracting the data for each location and putting it into 2D matrices
   [p1 p2 p3]=size(xxx);

% xx=zeros(p2,p3);
% xx1=zeros(p2,p3);
% yy=zeros(p2,p3);
% yy1=zeros(p2,p3); 

n=1;
n2=length(p);
%n2=9;
 
for k=n :n2

    k
    
    
    [p1 p2 p3]=size(xxx);
clear xx yy xx1 yy1
    
    xx=zeros(p2,p3);
xx1=zeros(p2,p3);
yy=zeros(p2,p3);
yy1=zeros(p2,p3); 


    for pp2=1:p2
        for pp3=1:p3
            xx(pp2,pp3)=xxx(k,pp2,pp3);
xx1(pp2,pp3)=xxx1(k,pp2,pp3);           
            yy(pp2,pp3)=yyy(k,pp2,pp3);
            yy1(pp2,pp3)=yyy1(k,pp2,pp3);
        end
    end
    
xx=transpose(xx);
xx1=transpose(xx1);

yy=transpose(yy);
yy1=transpose(yy1);

% close all
% plot(yy)
% figure
% plot(yy1)

nn=10;

for nnn = 1:nn         %start dialog to open file
    clear data_init data cp Fit1 dataf Nfft  fftm fftp Fs freq fftm_mask


% for object data 
x=xx1(:,nnn);
y=yy1(:,nnn)-yy(:,nnn) ;
    
   



% data operations
% 
% if k>7 && k<10
%     startind = find((y./max(y))>0.1,1)
%     
% else 
%    startind = find((y./max(y))>0.1,1)
% end
% if k>18 && k<21
%     startind = find((y./max(y))>0.15,1)
%     
% else 
%    startind = find((y./max(y))>0.1,1)
% end


   startind = find((y./max(y))>0.10,1)
% 
%    if k==75
% 
%    startind = find((y./max(y))>0.2,1)
% 
%    end
   
   % mar 
    [d ind] = max(y);
    data(:,1) = x(ind:end)-x(startind)+2.5*10^-10; %time choose T0
    
    % mar 
    %data(:,1) = x(ind:end)-x(startind)+6*10^-10; %time choose T0
    
    data(:,2) = y(ind:end); %amplitude
    dataplot(:,1) = x-x(startind)+2.5*10^-10; %time choose T0
    
    % mar 
    %dataplot(:,1) = x-x(startind)+6*10^-10; %time choose T0
    dataplot(:,2) = y; %amplitude
    
    
    if splot==1
    figure(1);
    hold off;
    plot(dataplot(startind-50:startind+500,1), dataplot(startind-50:startind+500,2));
    %mar
   % plot(dataplot(startind:startind+500,1), dataplot(startind:startind+500,2));
  %mar
    grid on
    hold on;
    plot(0, dataplot(startind,2), 'rd')
    grid on
    hold off;
    pause
    
    figure(2);
    hold off
    plot(dataplot(30:1800,1),dataplot(30:1800,2),'-b');
    grid on
   pause
    end
  

    
    %% do a rough THERMAL fit...
    
    LBp = [0 0 0]; %lower bound for fit variable in weird alphabetical order -
    %uppercase first [A B tau_th tau_zh time0] according to fit function below
    
    STp = [(max(data(:,2))) 13000 (mean(data(end-100:end,2)))]; %starting value for fit variable in
    %weird alphabetical order - uppercase first [A tau_th t0]
    %good starting values are needed for the fit to converge
    
    OPTIONS = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',LBp,...
        'Startpoint',STp);
    %             'Display','iter');  %these specifies type of fitting scheme and bounds
    %  and starting values for all fit variables
    
    TYPE=fittype('Ap.*erfc(Bp.*sqrt(x))+ Gp;','options',OPTIONS);
    [cp,gofp] = fit(data(:,1),data(:,2),TYPE);  %fit data - don't use ';' 
    %if you want to see iterations and final fit values
    
    %calculate the fit...
    Ap = cp.Ap;
    Bp = cp.Bp;
    Gp = cp.Gp;
    
    Fitrough = Ap.*erfc(Bp.*sqrt(data(:,1)))+Gp;
    
%        diffus_rough(nnn)=(Bp).^2./((2*pi/2.758e-6)^2);
       diffus_rough(nnn)=(Bp).^2./((2*pi/5.7289e-6)^2);
%     figure(2);
%     hold on
%     plot(data(:,1),Fitrough, '-r');
%     hold off
%     
%     pause
%     
%     

% taking the sinusiod part only - for the SAW fit 
    dataf = data(:,2)- Fitrough;
         
    
    %% take fft to determine sinusoid frequency...
    
    % makign lenght liek 1024 uif it was 1020 since fft is better  when
    % length is a multiple of 2
    Nfft = 2^(nextpow2(size(data(:,2),1)+1)+4); %size of zero padding for fft  , this gives smoother maximums 
      Y = fft(dataf, Nfft);   % matlab doesthe zero padding
     fftm = abs(Y((1:Nfft/2+1)));   % taking half of the fourirer nad getting the abs
    fftp = angle(Y((1:Nfft/2+1)));   % getting hte phase 
     
    Fs = 1/(data(end,1)-data(1,1))*(size(dataf,1)-1);   % makimg the frequency axis 
    freq = Fs/2*linspace(0,1,Nfft/2+1);
    
    fftm_mask = fftm; fftm_mask(1:150,:) = 0;     % removign the DC component max at zero ferequency 
    [max_int,fcI] = max(fftm_mask);               % now taking hte maximum  
    fc = freq(fcI); %centre frequency
    pc = fftp(fcI); %peak phase
    
    saw_f(nnn)=fc;
    %% thermal and SAW fit    
    %% Consider fits with different start positions... 
%     figure 
%     plot(freq,fftm_mask)
%     grid on 
%     xlabel('Freq. (Hz)')
%    
    %for iii = 1: ceil(1/fc/(data(2,1)-data(1,1))/1);     %  doing for the numnber wihtin one SAW period 
   for iii = 1:22    ;     %  doign for the numnber wihtin one SAW period
       % 12 is half saw 22 is full saw approxx
    
   
   iii
   
   delay=17;
   
  
        fitstart = (iii-1)*2 + delay  ;      % it does in jumps of 2                                
         
        %fitstart=40;
        fitend = 1700; %800 for 2.75 um grating %2500 for 9 um grating   
        % just the length of the trace we're fitting 
        
        % for the phase 
        
%         
%         if k>7 && k<10
%         
%         ph=-pi/3 -  0.1;
%         else  if k>18 && k<21
%             ph=-pi/3 - 0.1 ;
%             else
%                 ph=-pi/3;
%             end
%         end
%         
        ph=-pi/3+0.5;
      
        
%         ph=ph-(iii/22)*2*pi;
%         fc=1e12;
        
        %% Fit thermal profile:
        x0 = [(max(data(:,2))/erfc(1)/2) 13000 (max(data(:,2))/erfc(1)/2)/10 fc ph 10^-7 (mean(data(end-100:end,2)))]; %starting value for fit variable in weird alphabetical order - uppercase first [A tau_th t0]
        [yfit,x,~,~,~,~]=leasqr(data(fitstart:fitend,1),data(fitstart:fitend,2),x0,'decay_inc_ampl',0.0001,40);   
        
        
        %% x0 is the intiial gues for parameters, decay is the function , 
        %0.0001 is the tolerance and 40 is the iterations 
        
        %reconstruct the fitted data
        Ap_l(iii) = x(1);
        Bp_l(iii) = x(2);
        Cp_l(iii) = x(3);
        Dp_l(iii) = x(4);
        Ep_l(iii) = x(5);
        Fp_l(iii) = x(6);
        Gp_l(iii) = x(7);
        
        
        % at the end of each fit, saving hte parameters
        
        % plotting that fit  % some are one 1 some are L -l - change it 
        
        if iiiplot==1
         figure
         Fit2 = Ap_l(iii).*erfc(Bp_l(iii).*sqrt(data(fitstart:fitend,1)))+ Cp_l(iii).*cos((Dp_l(iii)*2*pi).*data(fitstart:fitend,1) + Ep_l(iii)).*exp(-data(fitstart:fitend,1)./Fp_l(iii)) + Gp_l(iii);
     plot(data(fitstart:fitend,1),Fit2, '-k');
     grid on
      hold on 
       plot(dataplot(30:1800,1),dataplot(30:1800,2),'-b');
        pause 
        end
    end
    
    %% taking the mean of the parameters 
    
    Ap_mean = mean(Ap_l);
    Bp_mean_l(nnn)= mean(Bp_l); 
    Bp_std_l(nnn)= std(Bp_l);
    Cp_mean = mean(Cp_l);
    fc_mean_l(nnn)= mean(Dp_l);
    fc;
    Ep_mean = mean(Ep_l);
    Fp_mean = mean(Fp_l);
    Gp_mean = mean(Gp_l);
    
    fitstart =(iii-1)*1 + 40; 
   % fitstart=40;
    
   clear Ap_l Bp_l Cp_l Dp_l Ep_l Fp_l Gp_l  
   
    % % using those to get the master fit 
    %full fit including sinusoid
    Fit1 = Ap_mean.*erfc(Bp_mean_l(nnn).*sqrt(data(fitstart:fitend,1)))+ Cp_mean.*cos((fc_mean_l(nnn)*2*pi).*data(fitstart:fitend,1) + Ep_mean).*exp(-data(fitstart:fitend,1)./Fp_mean) + Gp_mean;
     
    

%% plotting the master fit on the data     

if nnnplot==1
    figure(2);
    plot(dataplot(30:1800,1),dataplot(30:1800,2),'-b');
    hold on
    plot(data(fitstart:fitend,1),Fit1, '-k');
    grid on
% 
       pause
 end
    close all
    nnn
   
 end
% %% calculating the diffusivity and plotting
% 
% [Bp_mean_l', Bp_std_l'];

  diffus=(Bp_mean_l).^2./((2*pi/2.758e-6)^2);
 % diffus=(Bp_mean_l).^2./((2*pi/3.9481e-6)^2);

% diffus=(Bp_mean_l).^2./((2*pi/5.7289e-6)^2);
%  
%      fc_mean_l;
%     

%%
% %doing online calibration - tungsten only
%     mean(saw_f);
%     L=2670/mean(saw_f);
%     diffus=(Bp_mean_l).^2./((2*pi/L)^2);
    
%%    
  % time2=time(10:end);
%    plot(time2,diffus(1:end),'b')
%  
%        set(gcf,'color','w');
% set(gca,'fontsize',12);
% grid on
% xlabel('Time (s)','FontSize',14)
% ylabel('Thermal Diffusivity (Wm^-1K^-1)','FontSize',14)
   
%plot(time(1:n)/60,diffus)
% 
% plot(diffus)
% grid on
% xlabel('Test Time (min)')
% ylabel('Diffusivity Wm-1K-1')
% title('map1 ')
% 
% figure 
% %plot(time(1:n)/60,fc_mean_l)
% plot(fc_mean_l)
% 
% grid on
% xlabel('Test Time (min)')
% ylabel('SAW Freq')
% title('map1')

% 
map_diffuse(k)=mean(diffus);
std_diffuse(k)=std(diffus);
% map_diffuse(k)=mean(diffus_rough);
% std_diffuse(k)=std(diffus_rough);

% map_saw(k)=mean(fc_mean_l);
% 
% std_saw(k)=std(fc_mean_l);
map_saw(k)=mean(saw_f);
 std_saw(k)=std(saw_f);
 end    

%% plotting  
close all

figure 

% n=1;
% n2=85;
% p=p-0.5;
errorbar(p(n:n2),map_diffuse(n:n2),std_diffuse(n:n2))

hold on
plot(p(n:n2),map_diffuse(n:n2), 'r','LineWidth',1)
grid on 
xlabel('Location on Sample (mm)','FontSize',14)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',12);

% 
% saving=0;
if saving==1
savefig('sample3_line2_diffus_17delay.fig')
saveas(gcf,'sample3_line2_diffus_17delay.png')
end
figure 
% 
% 
errorbar(p(n:n2),map_saw(n:n2),std_saw(n:n2))

hold on
plot(p(n:n2),map_saw(n:n2), 'r','LineWidth',1)

grid on
xlabel('Location on Sample (mm)','FontSize',14)
ylabel('Peak SAW Frequency (Hz)','FontSize',14)
set(gcf,'color','w');
set(gca,'fontsize',12);
grid on

if saving==1
savefig('sample3_line2_SAW_17delay.fig')
saveas(gcf,'sample3_line2_SAW_17delay.png')

% 

% 
  save sample3_line2_analysis_17delay.mat p map_diffuse std_diffuse map_saw std_saw xxx xxx1 yyy yyy1

end

% 
% figure 
% 
% plot(map_diffuse(n:end), 'r','LineWidth',1)
% grid on 
% xlabel('Location on Sample (mm)','FontSize',14)
% ylabel('Thermal Diffusivity (Wm-1K-1)','FontSize',14)
% set(gcf,'color','w');
% set(gca,'fontsize',12);

%% cleaning plots 

% clear a
% z=6;
% 
% for i=1:10
%     for j=1:2002
%         a(i,j)=yyy1(z,i,j);
%     end
% end
% close all
% plot(a')
% grid on
