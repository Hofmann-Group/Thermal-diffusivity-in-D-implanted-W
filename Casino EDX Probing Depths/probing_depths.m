clc
clear all
close all

% reading the data file and taking the depth data and x-ray intesnity data
% (xray absorebed by detector)

% code loads all the data for the 5 MV case into variable W5MV
% then we extract depth of the 1/e of max intensity into W5MV1e, where 1e
% means the 1/e location
% MV denotes the MV characteristic x-ray, simillary at higher energies
% there is the L3 chracteristic xray
% the probing depths are then output into the command window, and the depth
% profiles are plot 

W5MV=dlmread('tungsten_5keV_MV_xray_distribution.dat','',2,0);

% getting the depth at 1/e of max intensity 
W5MV1e=W5MV((find(W5MV(:,3)./max(W5MV(:,3))<0.368,1)),1);



figure 
plot(W5MV(:,1),W5MV(:,3),'k^')
grid on 
xlabel('Depth (nm)','FontSize',16)
ylabel('Intensity ','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);

hold on


W10MV=dlmread('tungsten_10keV_MV_xray_distribution.dat','',2,0);
plot(W10MV(:,1),W10MV(:,3),'bd')


% getting the depth at 1/e of max intensity 
W10MV1e=W10MV((find(W10MV(:,3)./max(W10MV(:,3))<0.368,1)),1);


hold on

W15MV=dlmread('tungsten_15keV_MV_xray_distribution.dat','',2,0);
plot(W15MV(1:end-1,1),W15MV(1:end-1,3),'rx')



% getting the depth at 1/e of max intensity 
W15MV1e=W15MV((find(W15MV(:,3)./max(W15MV(:,3))<0.368,1)),1);


hold on

W15l3=dlmread('tungsten_15keV_l3_xray_distribution.dat','',2,0);
plot(W15l3(:,1),W15l3(:,3),'mo')

% getting the depth at 1/e of max intensity 
W15l31e=W15l3((find(W15l3(:,3)./max(W15l3(:,3))<0.368,1)),1);


hold on


W20MV=dlmread('tungsten_20keV_MV_xray_distribution.dat','',2,0);
plot(W20MV(:,1),W20MV(:,3),'c>')

% getting the depth at 1/e of max intensity 
W20MV1e=W20MV((find(W20MV(:,3)./max(W20MV(:,3))<0.368,1)),1);

hold on

W20l3=dlmread('tungsten_20keV_l3_xray_distribution.dat','',2,0);
plot(W20l3(:,1),W20l3(:,3),'g*')

% getting the depth at 1/e of max intensity 
W20l31e=W20l3((find(W20l3(:,3)./max(W20l3(:,3))<0.368,1)),1);



axis([0 500 0 3.5 ])
legend('5 keV MV','10 keV MV','15 keV MV','15 keV LIII','20 keV MV','20 keV LIII' )



depth=[W5MV1e  W10MV1e W15MV1e W15l31e W20MV1e W20l31e]
