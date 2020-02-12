% A basic test script demonstrating how to run the unifoil model.
% All units are in radians except when otherwise stated (eg some basis
% parameters)

%% Setup
clear;close all

%% Simulation Options
T = 500; % Simulation duration

%% Kite properties
% Physical properties
baseMass        = 6184;
addedMass       = 739.6;
baseInertia     = 80302;%104850;
addedInertia    = 0;%724530;
buoyFactor      = 1.05;
ARefWing        = 20;
ARefRudder      = 1.875;
ARefElevator    = 1.875;
fuselageLength  = 8;
wingOE          = 0.8;
rudderOE        = 0.8;
wingAR          = 5;
rudderAR        = 3;
wingTable       = buildAirfoilTable('wing',wingOE,wingAR);
rudderTable     = buildAirfoilTable('rudder',wingOE,wingAR);

% Initial Conditions
radius          = 100;
initSpeed       = 5.735;
initAzimuth     = 0.01*pi/180;
initElevation   = 29.99*pi/180;
initTwist       = -22.5*pi/180;
initTwistRate   = 0;

%% Water properties
flowSpeed   = 1;
density     = 1000;

%% Controller parameters
tauRef          = 0.025; % reference model time constant (s)
azimuthSweep    = 60; % Path azimuth sweep angle, degrees
elevationSweep  = 10; % Path elevation sweep angle, degrees
refGain1        = 1/tauRef^2;
refGain2        = 2/tauRef;
meanAzimuth     = 0;
meanElevation   = 30;
initPathVar     = 0;
maxLeadLength   = 0.01;
maxIntAngle     = 3*pi/180;
% Min and max alpha for the wing controller
wingAlphaPlusStall  = 10*pi/180;
wingAlphaMinusStall = -10*pi/180;
% Min and max alpha for the rudder controller
rudderAlphaPlusStall    = 6*pi/180;
rudderAlphaMinusStall   = -6*pi/180;
% Coefficients of linear and quadratic fits to CL and CD curves
[wingCLCoeffs,wingCDCoeffs]     = fitTable(wingTable,5*[-1 1]*pi/180);
[rudderCLCoeffs,rudderCDCoeffs] = fitTable(rudderTable,6*[-1 1]*pi/180);

%% Run the simulation
tic
sim('unifoil');
fprintf('Sim Efficiency: %.1f x Real Time\n',T/toc)

%% Post-process Data (get it into a signalcontainer object)
tsc = signalcontainer(logsout);

%% Find first time step of second lap
tIndx = find(and(tsc.pathVar.Data(1:end-1)>0.99,tsc.pathVar.Data(2:end)<0.01));
tVal  = tsc.pathVar.Time(tIndx(1));
tsc.crop([0 tVal]);

%% Plot some things
basisParams = [azimuthSweep, elevationSweep, meanAzimuth, meanElevation, radius];
path        = lemOfGerono(linspace(0,1),basisParams);
tsc.posVec.plot3('LineWidth',1,'Color','b','LineStyle','-','DisplayName','Flight Path')
daspect([1 1 1])
hold on
grid on
scatter3(0,0,0,'MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Origin')
plot3(path(:,1),path(:,2),path(:,3),...
    'LineWidth',2,'Color','r','LineStyle',':','DisplayName','Target Path')
view(84,40)
legend

figure
tsc.pathVar.plot