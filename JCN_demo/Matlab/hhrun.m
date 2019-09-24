function [V,m,n,h,INa,IK,Il]=hhrun(I,t)
%Constants set for all Methods
% dt=0.04; % Time Step ms
% t=0:dt:25; %Time Array ms
% I=0.1; %External Current Applied
dt=t(2)-t(1);

% Cm=0.01; % Membrane Capcitance uF/cm^2
% ENa=55.17; % mv Na reversal potential
% EK=-72.14; % mv K reversal potential
% El=-49.42; % mv Leakage reversal potential
% gbarNa=1.2; % mS/cm^2 Na conductance
% gbarK=0.36; % mS/cm^2 K conductance
% gbarl=0.003; % mS/cm^2 Leakage conductance
% V(1)=-60; % Initial Membrane voltage

% % params from Gerstner EPFL page
% Cm=1/200000; % uF/cm^2 / uF/mm2 divisé par 100, etc
% ENa=115-65; % mv Na reversal potential
% EK=-12-65; % mv K reversal potential
% El=10.6-65; % mv Leakage reversal potential
% gbarNa=120/200000; % mS/cm^2 Na conductance
% gbarK=36/200000; % mS/cm^2 K conductance
% gbarl=0.3/200000; % mS/cm^2 Leakage conductance
% V(1)=-65; % Initial Membrane voltage

% params from Gerstner EPFL page (http://icwww.epfl.ch/~gerstner/SPNM/node14.html)
Cm=1; % uF/cm^2 / uF/mm2 divisé par 100, etc
ENa=115-65; % mv Na reversal potential
EK=-12-65; % mv K reversal potential
El=10.7-65; % mv Leakage reversal potential
gbarNa=120; % mS/cm^2 Na conductance
gbarK=36; %36% mS/cm^2 K conductance
gbarl=0.3; % mS/cm^2 Leakage conductance
V(1)=-65; % Initial Membrane voltage

m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value
for i=1:length(t)-1
    %Euler method to find the next m/n/h value
    m(i+1)=m(i)+dt*((am(V(i))*(1-m(i)))-(bm(V(i))*m(i))); 
    n(i+1)=n(i)+dt*((an(V(i))*(1-n(i)))-(bn(V(i))*n(i)));
    h(i+1)=h(i)+dt*((ah(V(i))*(1-h(i)))-(bh(V(i))*h(i)));
    gNa=gbarNa*m(i)^3*h(i);
    gK=gbarK*n(i)^4;
    gl=gbarl;
    INa(i)=gNa*(V(i)-ENa);
    IK(i)=gK*(V(i)-EK);
    Il(i)=gl*(V(i)-El);
    %Euler method to find the next voltage value
    V(i+1)=V(i)+(dt)*((1/Cm)*(I(i)-(INa(i)+IK(i)+Il(i))));
end
INa(i+1)=gNa*(V(i+1)-ENa);
IK(i+1)=gK*(V(i+1)-EK);
Il(i+1)=gl*(V(i+1)-El);
end

% %Store variables for graphing later
% FE=V;
% FEm=m;
% FEn=n;
% FEh=h;
% clear V m n h;

% function a=am(v) 
% %Alpha for Variable m
% a=0.1*(v+35)/(1-exp(-(v+35)/10));
% end
% 
% function b=bm(v) 
% %Beta for variable m
% b=4.0*exp(-0.0556*(v+60));
% end
% 
% function a=an(v)
% %Alpha for variable n
% a=0.01*(v+50)/(1-exp(-(v+50)/10));
% end
% 
% function b=bn(v) 
% %Beta for variable n
% b=0.125*exp(-(v+60)/80);
% end
% 
% function a=ah(v) 
% %Alpha value for variable h
% a=0.07*exp(-0.05*(v+60));
% end
% 
% function b =bh(v) 
% %beta value for variable h
% b=1/(1+exp(-(0.1)*(v+30)));
% end

%% Gerstner page EPFL
function a=am(v) 
%Alpha for Variable m
v=v+65;
a=(2.5-0.1*v)/(exp(2.5-0.1*v)-1);
if v==25,
    a=1/2*((2.5-0.1*(v-1))/(exp(2.5-0.1*(v-1))-1)+(2.5-0.1*(v+1))/(exp(2.5-0.1*(v+1))-1));
end
end

function b=bm(v) 
%Beta for variable m
v=v+65;
b=4*exp(-v/18);
end

function a=an(v)
%Alpha for variable n
v=v+65;
a=(0.1-0.01*v)/(exp(1-0.1*v)-1);
if v==10,
    a=1/2*((0.1-0.01*(v-1))/(exp(1-0.1*(v-1))-1)+(0.1-0.01*(v+1))/(exp(1-0.1*(v+1))-1));
end
end

function b=bn(v) 
%Beta for variable n
v=v+65;
b=0.125*exp(-v/80);
end

function a=ah(v) 
%Alpha value for variable h
v=v+65;
a=0.07*exp(-v/20);
end

function b =bh(v) 
%beta value for variable h
v=v+65;
b=1/(exp(3-0.1*v)+1);
if v==30,
    b=1/2*(1/(exp(3-0.1*(v-1))+1)+1/(exp(3-0.1*(v+1))+1));
end
end
