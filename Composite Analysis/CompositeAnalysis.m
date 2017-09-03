function [ Sublaminate ] = CompositeAnalysis( FiberMaterial,MatrixMaterial,FvF,h,orient )
%COMPOSITEANALYSIS is a function that uses the defined Fiber and matrix
%materials as well as other input parameters to develop a laminate and
%analyze the Engineering Constants of it.
%
%INPUTS:
%   FiberMaterial = this is a cell array containing the fiber material
%       shorthand names
%   MatrixMaterial = this is a cell array containing the matrix material
%       shorthand names
%   FvF = this is a vector containing the Fiber Volume Fractions for each
%       lamina of the composite
%   h = this is a constant that defines the lamina thickness;
%   orient = this is a vector containing the fiber direction for each
%       lamina. Note that the angles are given in degrees.
%
%OUTPUTS:
%   Sublaminate = a structure containing important parameters of the
%       developed laminate. The fields of this structure are Ex, Ey, Gxy, 
%       vxy,vyx, ABD, Thickness, and Density;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numlayer=length(orient); % Determines the number of layers

for i=1:numlayer
    [ Fiber,Matrix ] = FindMaterial( FiberMaterial{i},MatrixMaterial{i} );
        % This function finds the fiber and matrix materials for each
        % layer. Both of these outputs are structures

    [ Composite(i) ] = CompositeProperties( Fiber,Matrix,FvF(i) );
        % This function uses the fiber and matrix materials and finds the
        % properties of the composite lamina. This output is a structure.
    Density(i)=Composite(i).Density*h;
        % This adds a field to the structure Composite that accounts for
        % the lamina density.
end

[Sublaminate]=Laminate(Composite,h,orient); % This function uses the Composite structure
    % as well as the layer thickness and the layer orientation vectors to 
    % develop the laminate and find the laminate properties
Sublaminate.Thickness=numlayer*h; % This finds the total thickness of the laminate
Sublaminate.Density=sum(Density)/Sublaminate.Thickness; % This finds the total density
    % of the laminate
end

function [ fiberprops,matrixprops ] = FindMaterial( fibermat,matrixmat )
FiberProperties;
MatrixProperties;

if isfield(Carbon,fibermat)
    fiberprops=eval(sprintf(['Carbon.',fibermat]));
elseif isfield(Glass,fibermat)
    fiberprops=eval(sprintf(['Glass.',fibermat]));
elseif isfield(Poly,fibermat)
    fiberprops=eval(sprintf(['Poly.',fibermat]));
else
    warning('Unrecognized material: %s',fibermat)
    custom=input('Would you like to enter custom fiber properties? [Y/N] :','s');
    switch custom
        case 'Y'
            fiberprops.YoungsModulus=input('Fiber Youngs Modulus = ');
            fiberprops.ShearModulus=input('Fiber Shear Modulus = ');
%             fiberprops.Poisson=input('Fiber Poissons Ratio = ');
        case 'N'
            error('Unrecognized material: \s',fibermat)
    end
end
if isempty(fiberprops.Poisson)
    fiberprops.Poisson=fiberprops.YoungsModulus/(2*fiberprops.ShearModulus)-1;
end
if isfield(Thermoplastic,matrixmat)
    matrixprops=eval(sprintf(['Thermoplastic.',matrixmat]));
elseif isfield(Thermoset,matrixmat)
    matrixprops=eval(sprintf(['Thermoset.',matrixmat]));
else
    warning('Unrecognized material: %s',matrixmat)
    custom=input('Would you like to enter custom matrix properties? [Y/N] :','s');
    switch custom
        case 'Y'
            matrixprops.YoungsModulus=input('Matrix Youngs Modulus = ');
            matrixprops.ShearModulus=input('Matrix Shear Modulus = ');
            matrixprops.Poisson=input('Matrix Poissons Ratio = ');
        case 'N'
            error('Unrecognized material: \s',matrixmat)
    end
end
if isempty(matrixprops.Poisson)
    matrixprops.Poisson=matrixprops.YoungsModulus/(2*matrixprops.ShearModulus)-1;
end
end

function [ Composite ] = CompositeProperties( Fiber,Matrix,FVF )
MVF=1-FVF;

Ef=Fiber.YoungsModulus;
Em=Matrix.YoungsModulus;
Gf=Fiber.ShearModulus;
Gm=Matrix.ShearModulus;
vf=Fiber.Poisson;
vm=Matrix.Poisson;

Composite.YoungsModulus.Principle1=FVF*Ef+MVF*Em;
Composite.YoungsModulus.Principle2=Ef*Em/(FVF*Em+MVF*Ef);    
Composite.YoungsModulus.Principle3=Composite.YoungsModulus.Principle2;

Composite.ShearModulus.Principle12=Gf*Gm/(FVF*Gm+MVF*Gf);   

Composite.Poisson.Principle12=vf*FVF+vm*MVF;

Composite.Density=FVF*Fiber.Density+MVF*Matrix.Density;
end

function [ laminate ] = Laminate( Composite,h,orient )
n=length(orient);
H=n*h;

N=length(orient);

z=zeros(1,N+1);
for i=1:N/2;
z(i)=-(N/2+1-i)*h;
z(N+2-i)=(N/2+1-i)*h;
end

[A,B,D,ABD]=CLT(z,orient,'deg',Composite);

a=A^(-1);

Ex=1/(a(1,1)*H);
Ey=1/(a(2,2)*H);
Gxy=1/(a(3,3)*H);
vxy=-a(1,2)/a(1,1);
vyx=-a(1,2)/a(2,2);

laminate.Ex=Ex;
laminate.Ey=Ey;
laminate.Gxy=Gxy;
laminate.vxy=vxy;
laminate.vyx=vyx;
laminate.ABD=ABD;
end

function [A,B,D,ABD]=CLT(z,Theta,cord,Composite)
for i=1:length(Theta)
    v=ones(1,3);
    G=ones(1,3);
    for j=1:3
        E(1,j)=getfield(Composite(i).YoungsModulus,sprintf(['Principle',num2str(j)]));
    end
    G(1,3)=Composite(i).ShearModulus.Principle12;
    v(1,3)=Composite(i).Poisson.Principle12;
    
    [S,C]=Compliance(E,v,G);
    [s,Q]=Reduce(S);
    
    Qs=[1 1 1/2;1 1 1/2;1 1 1/2];
    Ss=[1 1 2;1 1 2;1 1 2];
    
    T=Transformation(Theta(i),cord);
    Qt(:,:,i)=(T^(-1)*Q*T).*Qs;
    St(:,:,i)=(T^(-1)*s*T).*Ss;
end

A=zeros(size(Qt,1));
B=zeros(size(Qt,1));
D=zeros(size(Qt,1));
for k=1:(length(z)-1)
    A(:,:)=A(:,:)+(Qt(:,:,k).*(z(k+1)-z(k)));
    B(:,:)=B(:,:)+(1/2*Qt(:,:,k).*(z(k+1)^2-z(k)^2));
    D(:,:)=D(:,:)+(1/3*Qt(:,:,k).*(z(k+1)^3-z(k)^3));
end

ABD=[A B;B D];
end

function [ S,C ] = Compliance( e,V,G )
%INPUTS:
%   e=[E1 E2 E3]
%   V=[v23 v13 v12]
%   G=[G23 G13 G12]

E=[e(1) e(1) e(1) G(1) G(1) G(1);
   e(1) e(2) e(2) G(2) G(2) G(2);
   e(1) e(2) e(3) G(3) G(3) G(3);
   e(1) e(2) e(3) G(1) G(1) G(1);
   e(1) e(2) e(3) G(1) G(2) G(2);
   e(1) e(2) e(3) G(1) G(2) G(3)];

v=[  0   V(3) V(2)  0    0    0;
    V(3)  0   V(1)  0    0    0;
    V(2) V(1)  0    0    0    0;
     0    0    0    0    0    0;
     0    0    0    0    0    0;
     0    0    0    0    0    0];

I=eye(6);
S=zeros(6);
for i=1:6
    for j=1:6
        S(i,j)=S(i,j)+1/E(i,j)*I(i,j)-(1-I(i,j))*v(i,j)/E(i,j);
    end
end
C=S^(-1);
end

function [ s,Q ] = Reduce( S )
s=zeros(3);
Q=zeros(3);
q=zeros(3);

for i=1:2
    for j=1:2
        s(i,j)=S(i,j);
    end
end
s(3,3)=S(6,6)/2;

Q=s^(-1);
end

function out = Transformation( Theta,cord )
if strcmpi(cord,'deg')
    Theta=Theta*pi/180;
elseif strcmpi(cord,'rad')
    Theta=Theta;
else
    error('cord definition not valid')
end

out = [  cos(Theta)^2  sin(Theta)^2   2*cos(Theta)*sin(Theta)
         sin(Theta)^2  cos(Theta)^2  -2*cos(Theta)*sin(Theta)
        -cos(Theta)*sin(Theta)  cos(Theta)*sin(Theta)  cos(Theta)^2-sin(Theta)^2];
end


