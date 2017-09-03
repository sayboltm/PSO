
HEEDS=0; %This determines whether or not HEEDS is used as an optimizer (0=No,1=Yes)

LayerThickness=0.15e-3; % This determines the lamina thicknesses

% The following define the lamina Fiber Volume Fraction for each lamina
FvF1=0.5;
FvF2=0.5;
FvF3=0.5;
FvF4=0.5;
FvF5=0.5;
FvF6=0.5;

% Combines the FVFs for each lamina into a single vector
FvF=[FvF1 FvF2 FvF3 FvF4 FvF5 FvF6];

switch HEEDS
    case 0
        % The following define the lamina Fiber Material for each lamina
        FibMat1='PAN_UHM';
        FibMat2='PAN_UHM';
        FibMat3='Kevlar49';
        FibMat4='Kevlar49';
        FibMat5='PAN_UHM';
        FibMat6='PAN_UHM';
        
        % Combines the Fiber materials for each lamina into a single cell array
        FiberMaterial={FibMat1 FibMat2 FibMat3 FibMat4 FibMat5 FibMat6};
        
        % The following define the lamina Matrix Material for each lamina
        MtrxMat1='Polyester';
        MtrxMat2='Polyester';
        MtrxMat3='Epoxy';
        MtrxMat4='Epoxy';
        MtrxMat5='Polyester';
        MtrxMat6='Polyester';

        % Combines the Matrix materials for each lamina into a single cell array
        MatrixMaterial={MtrxMat1 MtrxMat2 MtrxMat3 MtrxMat4 MtrxMat5 MtrxMat6};

    case 1
        FibMat1=1;
        FibMat2=1;
        FibMat3=2;
        FibMat4=2;
        FibMat5=1;
        FibMat6=1;
        
        fibermaterial=[FibMat1 FibMat2 FibMat3 FibMat4 FibMat5 FibMat6];
        
        for i=1:length(fibermaterial)
            if fibermaterial(i)==1
                FiberMaterial{i}='PAN_UHM';
            elseif fibermaterial(i)==2
                FiberMaterial{i}='Kevlar49';
            elseif fibermaterial(i)==3
                FiberMaterial{i}='Graphite';
            elseif fibermaterial(i)==4
                FiberMaterial{i}='PAN_IM';
            elseif fibermaterial(i)==5
                FiberMaterial{i}='PAN_HM';
            elseif fibermaterial(i)==6
                FiberMaterial{i}='Pitch';
            elseif fibermaterial(i)==7
                FiberMaterial{i}='E';
            elseif fibermaterial(i)==8
                FiberMaterial{i}='Kevlar29';
            end
        end
        
        MtrxMat1=1;
        MtrxMat2=1;
        MtrxMat3=2;
        MtrxMat4=2;
        MtrxMat5=1;
        MtrxMat6=1;
        
        matrixmaterial=[MtrxMat1 MtrxMat2 MtrxMat3 MtrxMat4 MtrxMat5 MtrxMat6];
        
        for i=1:length(matrixmaterial)
            if matrixmaterial(i)==1
                MatrixMaterial{i}='Polyester';
            elseif matrixmaterial(i)==2
                MatrixMaterial{i}='Epoxy';
            elseif matrixmaterial(i)==3
                MatrixMaterial{i}='PEEK';
            end
        end
end

% The following define the number of sublaminates that are used in each
% section of the beam
Reg1=3;
Reg2=3;
Reg3=3;
Reg4=2;
Reg5=1;
Reg6=1;

% Combines the sublaminate numbers into a single vector
BeamRegions=[Reg1 Reg2 Reg3 Reg4 Reg5 Reg6];


