function [ Un,MaxStress,MaxDeflection,Mass,x ] = FiniteAnalysisExecution( Sublaminate,Regions )
%FINITEANALYSISEXECUTION Summary of this function goes here
%   Detailed explanation goes here

% if verLessThan('matlab','8.2.0')
%     if matlabpool('size')==0
%         set(findResource(),'ClusterSize',8);
%         matlabpool open 8
%     end
% else
%     if isempty(gcp('nocreate'))
%         parpool;
%     end
% end

[ x0,xL,BeamHeight,u0,uL,du0,duL,Mu0,MuL,Fu0,FuL,~ ] = InputFnc(); % Loads Beam information
[Z,Phi,dPhi,d2Phi]=Basis(); % This loads the basis functions

L=xL-x0; % Length of the beam
R=length(Regions); % Number of discrete regions of the beam. A region is a section of the beam where the cross section changes

%% Sublaminate Properties
Ex=Sublaminate.Ex; % Effective Young's Modulus of the composite sublaminate in the x direction
Ey=Sublaminate.Ey; % Effective Young's Modulus of the composite sublaminate in the x direction
vxy=Sublaminate.vxy; % Effective Major Poisson's ratio
vya=Sublaminate.vyx; % Effective Minor Poisson's ratio
p=Sublaminate.Density; % Density of the composite sublaminate
T=Sublaminate.Thickness;% Thickness of the sublaminate

%% Mesh
N=15*R; % Defines the number of elements to use in the FE analysis
x=linspace(x0,xL,N+1); % Creates a discretized domain of the beam.

c=zeros(N,4);
for i=1:N % This for loop develops the "basis connectivity" matrix (similar in concept to the nodal connectivity matrix)
    c(i,1)=2*i-1;
    c(i,2)=2*i;
    c(i,3)=2*i+1;
    c(i,4)=2*i+2;
end

%% Beam Regions
Mass=0; % Set the beam mass to 0
m=L/R; % finds the length of a region
j=1;
I=zeros(1,length(x));
Q=zeros(1,length(x));
for i=1:R % This loop finds the moments of inertias and other geometry based properties
    a=find(x>=i*m,1,'first'); % Find the first node in which the region changes
    for k=j:a % from the previous node to the first node in which the region changes
        b(k)=Regions(i)*T; % Determines the beam thickness for each element
        I(k)=1/12*b(k)*BeamHeight^3; % Finds the first (?) area moment of inertia for each element
        Q(k)=1/2*b(k)*((BeamHeight/2)^2); % Finds the second (?) area moment of inertia for each element
    end
    Mass=Mass+p*Regions(i)*T*BeamHeight*(x(a)-x(j)); % Determines the mass for the current region of the beam
    j=k+1;
end
EI=I.*Ex; % This is a vector containing the stiffness constant for each node

%% FE Analysis
parfor e=1:N % This loop finds the elemental entities 
    [z,dzdx]=Mapping(x,Z,e); % converts the element domain to isoparametric domain
    ei=EI(e); % sets the current stiffness constant
    Element(e).K=StiffnessMatrix(ei,dzdx,d2Phi); % finds the elemental stiffness matrix
    Element(e).R=ForcingVector(z,dzdx,Phi); % finds the elemental forcing vector
end

Le=x(2)-x(1);
K=8*EI(2)/Le^3.*[3/2 3/4*Le -3/2 3/4*Le; 0 1/2*Le^2 -3/4*Le 1/4*Le^2; 0 0 3/2 -3/4*Le; 0 0 0 1/2*Le^2];
Element(1).K;

Kg=sparse(2*N+2,2*N+2);
Rg=sparse(2*N+2,1);
for e=1:N % This for loop develops the global Stiffness and forcing vector from the elementary matrices and vectors
    Kg(c(e,:),c(e,:))=Kg(c(e,:),c(e,:))+Element(e).K;
    Rg(c(e,:),1)=Rg(c(e,:),1)+Element(e).R;
end

%% Boundary Condtions
Meq0=zeros(1,N);
Meqi=zeros(1,N);
MeqL=zeros(1,N);
Feq=zeros(1,N);
X=sym('x');
for e=1:N % Domain Force and Moment integration
    Mintegrand0=-(X-x(1))*fn(X);
    Mintegrandi=-(X-x(e))*fn(X);
    MintegrandL=-(x(e+1)-X)*fn(X);
    Feq(e)=Gaussian(fn(X),x(e),x(e+1),5);
    Meq0(e)=Gaussian(Mintegrand0,x(e),x(e+1),5);
    Meqi(e)=Gaussian(Mintegrandi,x(e),x(e+1),5);
    MeqL(e)=Gaussian(MintegrandL,x(e),x(e+1),5); 
end

P=[-1 1 0 0;0 -L -1 1]; % Boundary Condition Matrix
R=[sum(Feq) sum(Meq0)]'; % Boundary Condition Vector
resp=['R0';'RL';'M0';'ML']; % Boundary Condition element name

while (1) % Boundary Conditions
if ~isempty(u0) % If u0 is known...
    Kg(1,:)=0;
    Kg(1,1)=1;
    Rg(1,1)=u0;
end
if ~isempty(uL) % If uL is known...
    Kg(end-1,:)=0;
    Kg(end-1,end-1)=1;
    Rg(end-1,1)=uL;
end
if ~isempty(du0) % If du0 is known...
    Kg(2,:)=0;
    Kg(2,2)=1;
    Rg(2,1)=du0;
end
if ~isempty(duL) % If duL is known...
    Kg(end,:)=0;
    Kg(end,end)=1;
    Rg(end,1)=duL;
end
if ~isempty(Fu0) % If R1 is known...
    Rg(1)=Rg(1)+Fu0;
    R(:)=R(:)-P(:,end-3).*Fu0;
    P(:,end-3)=[];
    resp(end-3,:)=[];
    R0=Fu0;
end
if ~isempty(FuL)% If R2 is known...
    Rg(end-1)=Rg(end-1)+FuL;
    R(:)=R(:)-P(:,end-2).*FuL;
    P(:,end-2)=[];
    resp(end-2,:)=[];
    RL=FuL;
end
if ~isempty(Mu0) % If M1 is known...
    Rg(2)=Rg(2)+Mu0;
    R(:)=R(:)-P(:,end-1).*Mu0;
    P(:,end-1)=[];
    resp(end-1,:)=[];
    M0=Mu0;
end
if ~isempty(MuL)% If M2 is known...
    Rg(end)=Rg(end)+MuL;
    R(:)=R(:)-P(:,end).*MuL;
    P(:,end)=[];
    resp(end,:)=[];
    ML=MuL;
end
break
end

Response=P\R; % Finds the unkown responses (Forces and moments at boundaries)

for i=1:length(Response) % Defines the unknown responses
    if strcmpi(resp(i,:),'R0')
        R0=Response(i,:);
    elseif strcmpi(resp(i,:),'RL')
        RL=Response(i,:);
    elseif strcmpi(resp(i,:),'M0')
        M0=Response(i,:);
    else
        ML=Response(i,:);
    end
end

%% Response
Wn=Kg\Rg; % Finds the response of the beam to the applied load

Un=zeros(1,N+1);
Theta_n=zeros(1,N+1);

j=0;
k=0;
for i=1:length(Wn); % Develops the Displacement and rotation vectors from the beam response vector
    if mod(i,2)
        j=j+1;
        Un(j)=Wn(i); % This is the deflection vector
    else
        k=k+1;
        Theta_n(k)=Wn(i); % This is the slope vector
    end
end

%% Stresses and Strains
Strain=zeros(1,N);
F=zeros(1,length(x));
V=zeros(1,N+1);
M=zeros(1,N+1);
BendingStress=zeros(1,N+1);
ShearStress=zeros(1,N+1);
V(1)=R0;
M(1)=M0;
BendingStress(1)=M(1)*BeamHeight/(2*I(1));
ShearStress(1)=V(1)*Q(1)/(I(1)*b(1));
for e=1:N % Finds the strain, internal shear, internal moment, bending stress, and shear stress
    L0=x(e+1)-x(e); % The original length of each element
    l=sqrt((L0)^2+(Un(e+1)-Un(e))^2); % The deformed length of each element
    Strain(e)=abs(l-L0)/L0; % The axial strain of each element
    V(e+1)=V(e)+Feq(e);% The internal shear at each element
    M(e+1)=M(e)+(x(e+1)-x(e))*V(e+1)+Meqi(e); % The internal moment at each element
    BendingStress(e+1)=M(e+1)*BeamHeight/(2*I(e+1)); % The internal bending stress of each element
    ShearStress(e+1)=V(e+1)*Q(e+1)/(I(e+1)*b(e+1)); %  The internal shear stress of each element
end
TensileStress=Strain.*Ex;

MaxDeflection=max(abs(Un-u0));

MagBendingStress=sqrt(BendingStress.*BendingStress);
MagShearStress=sqrt(ShearStress.*ShearStress);

for i=1:length(x)
    F(i)=abs(fn(x(i)));
    MaxStress(i)=max(MagBendingStress(i),MagShearStress(i));
end

MaxStress=max(MaxStress);

Strain=Strain.*(10^6);
TensileStress=TensileStress.*(10^(-6));
BendingStress=BendingStress.*(10^(-9));
ShearStress=ShearStress.*(10^(-6));

un=Un.*(10^3);
F=F.*(10^(-3));
V=V.*(10^(-3));
M=M.*(10^(-3));

%%%% COMMENTING OUT FOR USING WITH PSO (and running bunch of times) %%%%%%%

%% Plotting and Command Window Display
% while (1)
%     ScreenPos=get(0,'MonitorPositions');
%     
%     figure;
%     subplot(3,1,1),plot(x,F,'r')
%     xlabel('x (meters)');
%     ylabel('Force (kN)');
%     title('Magnitude of Applied Load');
%     
%     subplot(3,1,2),plot(x,un,'b');
%     xlabel('x (meters)');
%     ylabel('U(x) (mm)');
%     Title=['Deflection Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     
%     subplot(3,1,3),plot(x,Theta_n,'b');
%     xlabel('x (meters)');
%     ylabel('\theta(x) (radians)');
%     Title=['Slope Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     set(gcf,'Position',[ScreenPos(1,1) ScreenPos(1,2) ScreenPos(1,3)/2 ScreenPos(1,end)]);
%     
%     Title=['Reaction Plots for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     pth=pwd;
%     plotfolder='\Plots\Figures\';
%     mkdir(pth,plotfolder)
%     file=fullfile(pth,plotfolder,Title);
%     saveas(gca,file);
%     
%     X=linspace(x0,xL,length(TensileStress));
%     
%     figure;
%     subplot(4,1,1),plot(X,TensileStress)
%     xlabel('x (meters)');
%     ylabel('Tensile Stress (MPa)');
%     Title=['Tensile Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     
%     subplot(4,1,3),plot(x,BendingStress)
%     xlabel('x (meters)');
%     ylabel('Bending Stress (GPa)');
%     Title=['Bending Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     
%     subplot(4,1,2),plot(x,ShearStress)
%     xlabel('x (meters)');
%     ylabel('Shear Stress (MPa)');
%     Title=['Shear Stress Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     
%     subplot(4,1,4),plot(X,Strain)
%     xlabel('x (meters)');
%     ylabel('Strain (\mum/mm)');
%     Title=['Strain Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     set(gcf,'Position',[ScreenPos(end,1)/2 ScreenPos(end,2) abs(ScreenPos(end,1))/2 ScreenPos(end,end)]);
%     
%     Title=['Resultant Plots for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     file=fullfile(pth,plotfolder,Title);
%     saveas(gca,file);
%     
%     figure
%     subplot(2,1,2),plot(x,M)
%     xlabel('x (meters)')
%     ylabel('Internal Moment (kNm)')
%     Title=['Internal Moment Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     
%     subplot(2,1,1),plot(x,V)
%     xlabel('x (meters)')
%     ylabel('Internal Shear (kN)')
%     Title=['Internal Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     title(Title);
%     grid on
%     set(gcf,'Position',[ScreenPos(end,1) ScreenPos(end,2) abs(ScreenPos(end,1))/2 ScreenPos(end,end)]);
%     
%     Title=['Moment-Shear Diagram for a Composite Beam with Laminate Regions of ',mat2str(Regions)];
%     file=fullfile(pth,plotfolder,Title);
%     saveas(gca,file);
%     break
% end % Develops several plots
% 
% fprintf('[R0 RL] = [%4.0f , %4.0f]\n',R0,RL)
% fprintf('[M0 ML] = [%4.0f , %4.0f]\n',M0,ML)
% maxdefindex=find(abs(Un-u0)==abs(MaxDeflection),1,'first');
% fprintf('Location of Maximum Deflection = %4.2f meters\n',x(maxdefindex))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function out=StiffnessMatrix(EI,dzdx,d2Phi)
out=zeros(length(d2Phi));
for i=1:length(d2Phi)
    for j=1:length(d2Phi)
        integrand=d2Phi(i)*EI*d2Phi(j)*dzdx^3;
        out(i,j)=Gaussian(integrand,-1,1,5);
    end
end
end

function out=ForcingVector(z,dzdx,Phi)
out=zeros(length(Phi),1);
for i=1:length(Phi)
    integrand=Phi(i)*fn(z)*(dzdx)^(-1);
    out(i,1)=Gaussian(integrand,-1,1,5);
end
end

function out=fn(x)
[ x0,xL,~,~,~,~,~,~,~,~,~,f ] = InputFnc();
L=xL-x0;

if isnumeric(x)
    out=-f/2*sin(2*pi/L*x-pi/2)-f/2;
%     out=-3e3*dirac(x-3);
else
    out=sym(sprintf([num2str(-f/2),'*sin(2*pi/',num2str(L),'*',char(x),'-pi/2)-',num2str(f/2)]));
%     out=sym(sprintf([num2str(-3e3),'*dirac(',char(x),'-3)']));
end
end

function [ z,dzdx] = Mapping( x,Z,e )
% [Z,Phi,dPhi,d2Phi]=Basis();
z=1/2*(1-Z)*x(e)+1/2*(1+Z)*x(e+1);
dzdx=2/(x(e+1)-x(e));
end

function [ Z,Phi,dPhi,d2Phi ] = Basis(  )
Z=sym('Z','real');
Phi=[1-3*((Z+1)./2).^2+2*((Z+1)./2).^3;
    (Z+1).*(1-(Z+1)./2).^2;
    3*((Z+1)./2).^2-2*((Z+1)./2).^3;
    (Z+1).*(((Z+1)./2).^2-((Z+1)./2))];
dPhi=diff(Phi);
d2Phi=diff(dPhi);
end