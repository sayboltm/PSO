function [ Un,MaxStress,MaxDeflection,Mass ] = FiniteAnalysisExecution( Sublaminate,Regions )
%FINITEANALYSISEXECUTION Uses the information stored the the strucutre
%"Sublaminate" as well as the vector "Regions" to perform the Finite
%Element Analysis of a stepped composite Beam modeled as an Euler Beam
%INPUTS:
%   Sublaminate = this is a structure containing information regarding the
%       sublaminate building blocks
%   Regions = this is a vector that specifies the number of sublaminates
%       used in the different regions of the beam
%OUTPUTS:
%   Un = is the nodal vector containing the displacement of the beam
%   MaxStress = this is the maximum stress in the beam
%   MaxDeflection = the maximum deflection experienced by the beam
%   Mass = this is the mass of the beam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verLessThan('matlab','8.2.0')
    if matlabpool('size')==0
        set(findResource(),'ClusterSize',8);
        matlabpool open 8
    end
else
    if isempty(gcp('nocreate'))
        parpool;
    end
end

[ x0,xL,BeamHeight,u0,uL,du0,duL,Mu0,MuL,Fu0,FuL,~ ] = InputFnc(); % Loads the domain as well as the Boundary Conditions

L=xL-x0;
R=length(Regions);

%% Sublaminate Properties
Ex=Sublaminate.Ex;
p=Sublaminate.Density;
T=Sublaminate.Thickness;

%% Mesh
N=10*L; %
x=linspace(x0,xL,N+1); %Discretizes the domain

c=zeros(N,4);
for i=1:N % This for-loop develops the nodal matrix
    c(i,1)=2*i-1;
    c(i,2)=2*i;
    c(i,3)=2*i+1;
    c(i,4)=2*i+2;
end

%% Beam Regions
Mass=0;
m=L/R;
j=1;
I=zeros(1,length(x));
Q=zeros(1,length(x));
for i=1:R % This for-loop calculates the thickness, moments of inertia, and the mass of the beam at each element
    a=find(x>=i*m,1,'first');
    for k=j:a
        b(k)=Regions(i)*T;
        I(k)=1/12*b(k)*BeamHeight^3;
        Q(k)=1/2*b(k)*((BeamHeight/2)^2);
    end
    Mass=Mass+p*Regions(i)*T*BeamHeight*(x(a)-x(j));
    j=k+1;
end
EI=I.*Ex; % This is a vector that is used as the material constant in the Stiffness Matrix calculations

%% FE Analysis
parfor e=1:N
    [z,dzdx,Phi,~,d2Phi]=Mapping(x,e); % Isoperametrices the beam elements
    ei=EI(e);
    Element(e).K=StiffnessMatrix(z,ei,dzdx,d2Phi); % Finds the elemental stiffness matrix
    Element(e).R=ForcingVector(z,dzdx,Phi); % Finds the elemental forcing vector
end

Kg=sparse(2*N+2,2*N+2);
Rg=sparse(2*N+2,1);
for e=1:N % This loop combines the elemental matrix and vector into global elements
    Kg(c(e,:),c(e,:))=Kg(c(e,:),c(e,:))+Element(e).K;
    Rg(c(e,:),1)=Rg(c(e,:),1)+Element(e).R;
end

%% Boundary Condtions
if ~isempty(u0)
    Kg(1,:)=0;
    Kg(1,1)=1;
    Rg(1,1)=u0;
end
if ~isempty(uL)
    Kg(end-1,:)=0;
    Kg(end-1,end-1)=1;
    Rg(end-1,1)=uL;
end
if ~isempty(du0)
    Kg(2,:)=0;
    Kg(2,2)=1;
    Rg(2,1)=du0;
end
if ~isempty(duL)
    Kg(end,:)=0;
    Kg(end,end)=1;
    Rg(end,1)=duL;
end
if ~isempty(Mu0)
    Rg(2)=Rg(2)+Mu0;
end
if ~isempty(MuL)
    Rg(end)=Rg(end)+MuL;
end
if ~isempty(Fu0)
    Rg(1)=Rg(1)+Fu0;
end
if ~isempty(FuL)
    Rg(end-1)=Rg(end-1)+FuL;
end
%% Response
Wn=Kg\Rg;

Un=zeros(1,N+1);
Theta_n=zeros(1,N+1);

for i=1:length(Wn)/2
    Un(i)=Wn(2*i-1);
    Theta_n(i)=Wn(2*i);
end

%% Stresses and Strains
Strain=zeros(1,N);
F=zeros(1,length(x));
M=zeros(1,N);
BendingStress=zeros(1,N);
ShearStress=zeros(1,N);
for e=1:N
    L0=x(e+1)-x(e);
    l=sqrt((L0)^2+(Un(e+1)-Un(e))^2);
    Strain(e)=abs(l-L0)/L0;
end
TensileStress=Strain.*Ex;
for e=1:N
    [z,~,~,~,~]=Mapping(x,e);
    Feq(e)=Gaussian(fn(z),z,5);
    Meq(e)=Gaussian(z*fn(z),z,5);
end
x_bar=sum(Meq)/sum(Feq);
reaction=[0 xL;1 1]^(-1)*[Feq*x_bar; Feq];
V=Feq;
V(1)=V(1)+reaction(1);

for e=1:N
    [z,~,~,~,~]=Mapping(x,e);
    M(e)=Gaussian(int(fn(z)+reaction(1)),z,5);
    BendingStress(e)=M(e)*BeamHeight/(2*I(e));
    ShearStress(e)=V(e)*Q(e)/(I(e)*b(e));
end

MaxDeflection=max(abs(Un))-min(abs(Un));
StressMax=[max(abs(TensileStress)),max(abs(BendingStress)),max(abs(ShearStress))];
MaxStress=max(StressMax);

Strain=Strain.*(10^6);
TensileStress=TensileStress.*(10^(-6));
BendingStress=BendingStress.*(10^(-9));
ShearStress=ShearStress.*(10^(-6));

for i=1:length(x)
    F(i)=abs(fn(x(i)));
end

un=Un.*(10^3);
F=F.*(10^(-3));

%% Plotting
figure;
subplot(3,1,1),plot(x,F,'r')
xlabel('x (meters)');
ylabel('Force (kN)');
title('Magnitude of Applied Load');

subplot(3,1,2),plot(x,un,'b');
xlabel('x (meters)');
ylabel('U(x) (mm)');
Title=['Deflection Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);

subplot(3,1,3),plot(x,Theta_n,'b');
xlabel('x (meters)');
ylabel('\theta(x) (radians)');
Title=['Slope Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);
grid on

Title=['Reaction Plots for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
pth=pwd;
plotfolder='\Plots\Figures\';
file=fullfile(pth,plotfolder,Title);
saveas(gca,file);

X=linspace(x0,xL,length(TensileStress));

figure;
subplot(4,1,1),plot(X,TensileStress)
xlabel('x (meters)');
ylabel('Tensile Stress (MPa)');
Title=['Tensile Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);
grid on

subplot(4,1,3),plot(X,BendingStress)
xlabel('x (meters)');
ylabel('Bending Stress (GPa)');
Title=['Bending Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);
grid on

subplot(4,1,2),plot(X,ShearStress)
xlabel('x (meters)');
ylabel('Shear Stress (MPa)');
Title=['Shear Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);
grid on

subplot(4,1,4),plot(X,Strain)
xlabel('x (meters)');
ylabel('Strain (\mum/mm)');
Title=['Strain Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
title(Title);
grid on

Title=['Resultant Plots for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
file=fullfile(pth,plotfolder,Title);
saveas(gca,file);

% maxdefindex=find(abs(Un)==MaxDeflection);
% fprintf('Location of Maximum Deflection = %4.2f meters\n',x(maxdefindex))



end

function [ x0,xL,BeamHeight,u0,uL,du0,duL,Mu0,MuL,Fu0,FuL,f ] = InputFnc()

x0=0;
xL=5;

BeamHeight=500e-3;

u0=0;
uL=0;

du0=[];
duL=[];

Mu0=0;
MuL=0;

Fu0=[];
FuL=[];

f=1e3;
end

function out=StiffnessMatrix(z,EI,dzdx,d2Phi)
out=zeros(length(d2Phi));
for i=1:length(d2Phi)
    for j=1:length(d2Phi)
        integrand=d2Phi(i)*EI*d2Phi(j)*dzdx^3;
        out(i,j)=Gaussian(integrand,z,5);
    end
end
end

function out=ForcingVector(z,dzdx,Phi)
out=zeros(length(Phi),1);
for i=1:length(Phi)
    integrand=Phi(i)*fn(z)*(dzdx)^(-1);
    out(i,1)=Gaussian(integrand,z,5);
end
end

function out=fn(x)
[ x0,xL,~,~,~,~,~,~,~,~,~,f ] = InputFnc();
L=xL-x0;
if isnumeric(x)
    out=-f/2*sin(2*pi/L*x-pi/2)-f/2;
else
    out=sym(sprintf([num2str(-f/2),'*sin(2*pi/',num2str(L),'*',char(x),'-pi/2)-',num2str(f/2)]));
end
end

function [ z,dzdx,Phi,dPhi,d2Phi ] = Mapping( x,e )
[Z,Phi,dPhi,d2Phi]=Basis();
z=1/2*(1-Z)*x(e)+1/2*(1+Z)*x(e+1);
dzdx=2/(x(e+1)-x(e));
end

function [ Z,Phi,dPhi,d2Phi ] = Basis(  )
Z=sym('Z','real');
Phi=[1-3*((Z+1)/2)^2+2*((Z+1)/2)^3;
    (Z+1)*(1-(Z+1)/2)^2;
    3*((Z+1)/2)^2-2*((Z+1)/2)^3;
    (Z+1)*(((Z+1)/2)^2-((Z+1)/2))];
dPhi=diff(Phi);
d2Phi=diff(dPhi);
end