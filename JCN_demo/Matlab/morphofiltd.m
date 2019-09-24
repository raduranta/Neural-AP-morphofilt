function w = morphofiltd(re,ordre,r0,r1,rN,rD,Cs)
%MORPHOFILT Summary of this function goes here
%   inputs: elpos = electrode positions (M x 3)
%           ordre = filter length (N+1 = nb of compartments on the axon + soma)
%           r0 = soma position (1 x 3) (center)
%           r1 = axon hillock position (1 x 3) (begining)
%           rN = last axon compartment position (1 x 3) (begining)
%           rD = tip of the equivalent dendrite
%           Cs = amplitude of the somatic dipole

%dk=rN-r0/ordre;
if nargin<5,
    error('Not enough arguments');
elseif nargin<6
    rD=r1;
    Cs=1;
elseif nargin<7
    Cs=1;
end
cond=0.33;
M=size(re,1);
rk=[linspace(r1(1),rN(1),ordre-1);linspace(r1(2),rN(2),ordre-1);linspace(r1(3),rN(3),ordre-1)]';
w=zeros(M,ordre);
for iel=1:M
    for ik=2:ordre-1
        w(iel,ik)=-(re(iel,:)-rk(ik-1,:))*(rk(ik,:)-rk(ik-1,:))'/(4*pi*cond*norm(re(iel,:)-rk(ik-1,:))^3);
    end
    w(iel,ordre)=-(re(iel,:)-rk(ordre-1,:))*(rk(ordre-1,:)-rk(ordre-2,:))'/(4*pi*cond*norm(re(iel,:)-rk(ordre-1,:))^3);
    w(iel,1)=-Cs*(re(iel,:)-r0)*(rD-r0)'/(4*pi*cond*norm(re(iel,:)-r0)^3);
end
end

