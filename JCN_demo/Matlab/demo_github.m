clear all
close all

%% HH
inmvm=3000; % index max on Vm in LFPy (3000 for synchronisation)
lVLFPy=8000;% signal length in LFPy
dt=10^(-3); % in ms
Nt=2^15;
D=Nt*dt;
t=[dt:dt:D]-dt;
n=length(t);

fe=1/dt;
f=0:fe/Nt:fe/2-fe/Nt;

I=(heaviside(t-1)-heaviside(t-31))*0.044/(2*pi*12.5*25)*10^8*10^-3;%*5.093;%0.15/(pi*12.5*12.5*2+2*pi*12.5*25)*10^8;
icur=1;
[Vm,m,n,h,INa,IK,Il]=hhrun(I,t); % pot membrane, proportionnel au courant des canaux ioniques (http://www.bem.fi/book/03/03.htm, 3.14)
Im=(INa+IK+Il)*(2*pi*12.5*25)/10^8*10^3;
[MVm,inMVm]=max(Vm);
%% BS neuron morphology

SL=25; % soma length (cylinder with the same diameter)

LA=1000; %axon length
DA=2; % %axon diameter

LD=50; %dendrite length 
DD=2; %dendrite diameter
phi=pi/2; % angle avec Oz
theta=pi; % angle with Ox (phi=pi/2,theta=pi) indicates opposite to the axon

%% load LFPy simulation result

Vlfpy=dlmread(['../Python/Vlfpy_BS_LA',num2str(LA),'_DA',num2str(DA),'_LD',num2str(LD),'_DD',num2str(DD),'demo.txt']);
Vmlfpy=dlmread(['../Python/Vm_BS_LA',num2str(LA),'_DA',num2str(DA),'_LD',num2str(LD),'_DD',num2str(DD),'demo.txt']);
Imlfpy=dlmread(['../Python/Im_BS_LA',num2str(LA),'_DA',num2str(DA),'_LD',num2str(LD),'_DD',num2str(DD),'demo.txt']);

%% figure check
figure
subplot(2,1,1)
plot(Vm([inMVm-inmvm+1:inMVm-inmvm+lVLFPy])) 
hold on
plot(Vmlfpy)

subplot(2,1,2)
plot(Im([inMVm-inmvm+1:inMVm-inmvm+lVLFPy]))
hold on
plot(Imlfpy)


%% filter parameters
dk=10; % axonal spatial sampling (~ nb of segments)
ordre=LA/dk+1;
r0=[0 0 0]; % soma position
r1=[SL/2 0 0]; % axon start position
rN=[SL/2+LA-dk 0 0]; % axon stop position (start of the last segment)
rd=norm(r1-r0)*[sin(phi)*cos(theta) sin(phi)*sin(theta) cos(phi)]; % dendrite end position, normalized
Cs=2; % somatic equivalent dipole amplitude
taus=23; % subsampling of the membrane current dk/taus = speed v)

%% electrodes
X=[-250:125:1250]';
Y=[250:-50:50]';
Z=0;

[eplosy,elposx,elposz]=meshgrid(Y,X,Z);
elpos=[elposx(:),eplosy(:),elposz(:)]

%% simulation
w = morphofiltd(elpos,ordre,r0,r1,rN,rd,Cs);
wup=upsample(w',taus)';

Vel=zeros(size(w,1),length(Im));
for iel=1:size(w,1),
    Vel(iel,:)=conv(Im,wup(iel,:),'same');
end
% cut
intervVm=[inMVm-inmvm-fix(size(wup,2)/2)+1:inMVm-inmvm-fix(size(wup,2)/2)+lVLFPy];
Vel2=Vel(:,intervVm);
% normalize
elsync=56;
Vel2=Vel2/norm(Vel2(elsync,:))*norm(Vlfpy(:,elsync));
%% plot grid 

cc=zeros(1,size(elpos,1));
t=dt:dt:dt*size(Vel2,2);
figure
cmap=colormap;
for ifil=1:size(elpos,1),
    subplot(5,13,ifil);
    plot(t,Vel2(ifil,:)-Vel2(ifil,1),'LineWidth',2)
    hold on
    plot(t,Vlfpy(:,ifil)-Vlfpy(1,ifil),'LineWidth',2)
    cc(ifil)=corr(Vel2(ifil,:)',Vlfpy(:,ifil));
    scatter(4,-2*10^-3,100,cmap(1+fix(size(cmap,1)*cc(ifil)),:),'filled')
    ylim([-5 5]*10^-3)
    if ifil>52
        text(2,-12*10^-3,[num2str((ifil-53)*125-250),'\mu','m'])
    end
    if rem(ifil,13)==1,
       text(-8,0*10^-3,[num2str(-fix(ifil/13)*50+250),'\mu','m'])
    end 
    axis off
end
colorbar('Position',[0.93 0.3 0.007 0.6],'FontSize',14)

fprintf('\n Mean correlation = %1.2f \n Min correlation = %1.2f  \n Max correlation = %1.2f \n',mean(cc),min(cc),max(cc))
