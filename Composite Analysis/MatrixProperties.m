
Thermoset.Polyester.Density=(1100+1500)/2;
Thermoset.Polyester.YoungsModulus=((1.2+4.5)/2).*(10^9);
Thermoset.Polyester.ShearModulus=((0.7+2)/2).*(10^9);
Thermoset.Polyester.Poisson=0.33;
% Thermoset.Polyester.TensileStrength=[40 90].*(10^6);
% Thermoset.Polyester.CompressiveStrength=[90 250].*(10^6);
% Thermoset.Polyester.Elongation=[2 5].*0.01;
% Thermoset.Polyester.ThermalExpansion=[60 200].*(10^(-6));
% Thermoset.Polyester.Conductivity=0.2;
% Thermoset.Polyester.SpecificHeat=[];

Thermoset.VinylEster.Density=1150;
Thermoset.VinylEster.YoungsModulus=((3+4)/2).*(10^9);
Thermoset.VinylEster.ShearModulus=[]*10^9;
Thermoset.VinylEster.Poisson=[];
% Thermoset.VinylEster.TensileStrength=[65 90].*(10^6);
% Thermoset.VinylEster.CompressiveStrength=127*10^6;
% Thermoset.VinylEster.Elongation=[1 5].*0.01;
% Thermoset.VinylEster.ThermalExpansion=53*10^(-6);
% Thermoset.VinylEster.Conductivity=[];
% Thermoset.VinylEster.SpecificHeat=[];

Thermoset.Epoxy.Density=(1100+1400)/2;
% Thermoset.Epoxy.YoungsModulus=((2+6)/2).*(10^9);
% Thermoset.Epoxy.ShearModulus=((1.1+2.2)/2).*(10^9);
Thermoset.Epoxy.YoungsModulus=6*(10^9);
Thermoset.Epoxy.ShearModulus=2.2*(10^9);
Thermoset.Epoxy.Poisson=0.33;
% Thermoset.Epoxy.TensileStrength=[35 130].*(10^6);
% Thermoset.Epoxy.CompressiveStrength=[100 200].*(10^6);
% Thermoset.Epoxy.Elongation=[1 8.5].*0.01;
% Thermoset.Epoxy.ThermalExpansion=[45 70].*(10^(-6));
% Thermoset.Epoxy.Conductivity=[0.1 0.2];
% Thermoset.Epoxy.SpecificHeat=[1250 1800];

Thermoset.Bismaleimide.Density=1320;
Thermoset.Bismaleimide.YoungsModulus=3.6*10^9;
Thermoset.Bismaleimide.ShearModulus=1.8*10^9;
Thermoset.Bismaleimide.Poisson=[];
% Thermoset.Bismaleimide.TensileStrength=[48 78].*(10^6);
% Thermoset.Bismaleimide.CompressiveStrength=200*10^6;
% Thermoset.Bismaleimide.Elongation=[1 6.6].*0.01;
% Thermoset.Bismaleimide.ThermalExpansion=49*10^(-6);
% Thermoset.Bismaleimide.Conductivity=[];
% Thermoset.Bismaleimide.SpecificHeat=[];

Thermoset.Polyimide.Density=(1430+1890)/2;
Thermoset.Polyimide.YoungsModulus=((3.1+4.9)/2).*(10^9);
Thermoset.Polyimide.ShearModulus=[]*10^9;
Thermoset.Polyimide.Poisson=[];
% Thermoset.Polyimide.TensileStrength=[70 120].*(10^6);
% Thermoset.Polyimide.CompressiveStrength=[]*10^6;
% Thermoset.Polyimide.Elongation=[1.5 3].*0.01;
% Thermoset.Polyimide.ThermalExpansion=90*10^(-6);
% Thermoset.Polyimide.Conductivity=[];
% Thermoset.Polyimide.SpecificHeat=[];

Thermoplastic.Ultem.Density=1270;
Thermoplastic.Ultem.YoungsModulus=3*10^9;
Thermoplastic.Ultem.ShearModulus=3*10^9;
Thermoplastic.Ultem.Poisson=[];
% Thermoplastic.Ultem.TensileStrength=105*10^6;
% Thermoplastic.Ultem.CompressiveStrength=140*10^6;
% Thermoplastic.Ultem.Elongation=60*0.01;
% Thermoplastic.Ultem.ThermalExpansion=62*10^(-6);
% Thermoplastic.Ultem.Conductivity=[];

Thermoplastic.Torlon.Density=1400;
Thermoplastic.Torlon.YoungsModulus=5*10^9;
Thermoplastic.Torlon.ShearModulus=5*10^9;
Thermoplastic.Torlon.Poisson=[];
% Thermoplastic.Torlon.TensileStrength=[95 185].*(10^6);
% Thermoplastic.Torlon.CompressiveStrength=276*10^6;
% Thermoplastic.Torlon.Elongation=[12 18].*0.01;
% Thermoplastic.Torlon.ThermalExpansion=36*10^(-6);
% Thermoplastic.Torlon.Conductivity=[];

Thermoplastic.PPS.Density=1340;
Thermoplastic.PPS.YoungsModulus=3.3*10^9;
Thermoplastic.PPS.ShearModulus=3.3*10^9;
Thermoplastic.PPS.Poisson=[];
% Thermoplastic.PPS.TensileStrength=[70 75].*(10^6);
% Thermoplastic.PPS.CompressiveStrength=110*10^6;
% Thermoplastic.PPS.Elongation=3*0.01;
% Thermoplastic.PPS.ThermalExpansion=[54 100].*(10^(-6));
% Thermoplastic.PPS.Conductivity=[];

Thermoplastic.PEEK.Density=1320;
Thermoplastic.PEEK.YoungsModulus=3.6*10^9;
Thermoplastic.PEEK.ShearModulus=1.4*10^9;
Thermoplastic.PEEK.Poisson=0.37;
% Thermoplastic.PEEK.TensileStrength=[92 100].*(10^6);
% Thermoplastic.PEEK.CompressiveStrength=[]*10^6;
% Thermoplastic.PEEK.Elongation=150*0.01;
% Thermoplastic.PEEK.ThermalExpansion=[]*10^(-6);
% Thermoplastic.PEEK.Conductivity=[];

Thermoplastic.PS.Density=1240;
Thermoplastic.PS.YoungsModulus=2.5*10^9;
Thermoplastic.PS.ShearModulus=2.5*10^9;
Thermoplastic.PS.Poisson=[];
% Thermoplastic.PS.TensileStrength=[70 75].*(10^6);
% Thermoplastic.PS.CompressiveStrength=[]*10^6;
% Thermoplastic.PS.Elongation=[50 100].*0.01;
% Thermoplastic.PS.ThermalExpansion=[56 100].*(10^(-6));
% Thermoplastic.PS.Conductivity=[];

Thermoplastic.PP.Density=900;
Thermoplastic.PP.YoungsModulus=((1+1.4)/2).*(10^9);
Thermoplastic.PP.ShearModulus=((1+1.4)/2).*(10^9);
Thermoplastic.PP.Poisson=[];
% Thermoplastic.PP.TensileStrength=[25 38].*(10^6);
% Thermoplastic.PP.CompressiveStrength=[]*10^6;
% Thermoplastic.PP.Elongation=300*0.01;
% Thermoplastic.PP.ThermalExpansion=110*10^(-6);
% Thermoplastic.PP.Conductivity=0.2;

Thermoplastic.Nylon.Density=1140;
Thermoplastic.Nylon.YoungsModulus=((1.4+2.8)/2).*(10^9);
Thermoplastic.Nylon.ShearModulus=((1.4+2.8)/2).*(10^9);
Thermoplastic.Nylon.Poisson=[];
% Thermoplastic.Nylon.TensileStrength=[60 75].*(10^6);
% Thermoplastic.Nylon.CompressiveStrength=34*10^6;
% Thermoplastic.Nylon.Elongation=[40 80].*0.01;
% Thermoplastic.Nylon.ThermalExpansion=90*10^(-6);
% Thermoplastic.Nylon.Conductivity=0.2;

Thermoplastic.PC.Density=((1060+1200)/2);
Thermoplastic.PC.YoungsModulus=((2.2+2.4)/2).*(10^9);
Thermoplastic.PC.ShearModulus=((2.2+2.4)/2).*(10^9);
Thermoplastic.PC.Poisson=[];
% Thermoplastic.PC.TensileStrength=[45 70].*(10^6);
% Thermoplastic.PC.CompressiveStrength=86*10^6;
% Thermoplastic.PC.Elongation=[50 100].*0.01;
% Thermoplastic.PC.ThermalExpansion=70*10^(-6);
% Thermoplastic.PC.Conductivity=0.2;
