% A basic test script demonstrating how to run the unifoil model.
% All units are in radians except when otherwise stated (eg some basis
% parameters)

%% Setup
clear;close all

%% Simulation Options
T = 100; % Simulation duration

%% Kite properties
% Physical properties
baseMass        = 6184;
addedMass       = 739.6;
baseInertia     = 10000;%80302;%104850;
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
initElevation   = 30*pi/180;
initTwist       = -20*pi/180;%-22.5*pi/180;
initTwistRate   = 0;

%% Water properties
flowSpeed   = 1;
density     = 1000;

%% Controller parameters
% Path geometry
azimuthSweep    = 60*pi/180; % Path azimuth sweep angle, degrees
elevationSweep  = 10*pi/180; % Path elevation sweep angle, degrees
meanAzimuth     = 0*pi/180;
meanElevation   = 30*pi/180;
pathShape  = 2; % 1 for ellipse, 2 for fig 8
basisParams = [azimuthSweep, elevationSweep, meanAzimuth, meanElevation, radius, pathShape];
% Reference model
tauRef          = 0.025; % reference model time constant (s)
% Model Ref Ctrl Gains
refGain1        = 1/tauRef^2;
refGain2        = 2/tauRef;
% Pure Pursuit Controller
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

%% Flexible Time ILC Parameters
pathStep = 0.005;
% Linear force parmeterization
FxParams = [1000 8];

% Waypoint Specification
posWayPtPathVars = [0.25 0.5 0.75 1];
wyptAzimuthTol   = 0.5*pi/180;
wyptElevationTol = 0.5*pi/180;
ctrlInputTol     = 0.5*pi/180;

% Calculate derived quantities
wyPtPosVecs     = lemOfGerono(posWayPtPathVars,basisParams);
wyPtAzimuth     = atan2(wyPtPosVecs(:,2),wyPtPosVecs(:,1));
wyPtElevation   = pi/2-acos(wyPtPosVecs(:,3)./sqrt(sum(wyPtPosVecs.^2,2)));
TPhi    = wyptAzimuthTol*ones(size(wyPtAzimuth));
TTheta  = wyptElevationTol*ones(size(wyPtElevation));

% % Plot the modeled Fx, sanity check
% twist = linspace(-180,180)*pi/180;
% az    = linspace(-80,80)*pi/180;
% [twist,az] = meshgrid(twist,az);
% Fx = ca.*cos(twist).^2.*exp(-(az./cb).^2);
% surf(twist*180/pi,az*180/pi,Fx)
% xlabel('Twist');ylabel('Azimuth');zlabel('Fx')

%% Run the simulation
tic
sim('unifoil');
fprintf('Sim Efficiency: %.1f x Real Time\n',T/toc)

%% Post-process Data (get it into a signalcontainer object)
tsc = signalcontainer(logsout);
tscTime = tsc;

%% Do operations for flexible time ILC

% Crop to just one lap
tIndx = find(and(tscTime.pathVar.Data(1:end-1)>0.99,tscTime.pathVar.Data(2:end)<0.01));
tVal  = tscTime.pathVar.Time(tIndx(1));
tscTime.crop([0 tVal]);

% Reparameterize signals in terms of the path variable
tscPath = reparameterize(tscTime);

% Resample all signals to different path steps
tscPath = tscPath.resample(0:pathStep:1);

% Create continuous, linear, path parmeterized model
[Acp,Bcp] = pathLinearize(tsc,basisParams,FxParams,tauRef,baseMass+addedMass);

% Discretize the continuous, linear, path-domain model
[Adp,Bdp] = disc(Acp,Bcp,pathStep);

% Lift the discrete, linear, path-domain model
[F,G] = lift(Adp,Bdp);

% Build the weighting matrix for the performance index
J = ones(1,numel(tscPath.rudderCmd.Time))*stateSelectionMatrix(3,5,numel(Adp.Time))*G;
Qphi   = wyptSelectionMatrix(1,5,...
    cnvrtPathVar2Indx(posWayPtPathVars,tscPath.pathVar.Data),...
    numel(tscPath.pathVar.Data));
Qtheta = wyptSelectionMatrix(2,5,...
    cnvrtPathVar2Indx(posWayPtPathVars,tscPath.pathVar.Data),...
    numel(tscPath.pathVar.Data));

% Setup large A and b matrices that we need for linprog
x0      = tscPath.stateVec.getdatasamples(1);
ALP     = [Qphi*G; -Qphi*G; Qtheta*G; -Qtheta*G];
b0LP    = [TPhi;TPhi;TTheta;TTheta];
b1LP    = [-Qphi*F*x0(:); Qphi*F*x0(:) ; -Qtheta*F*x0(:); Qtheta*F*x0(:) ];
bLP     = b0LP - b1LP;

% Run linprog to get the update to the control input
duStar = linprog(J,ALP,bLP);

% ILC Update
uNext = uPrev + duStar; % uPrev should be the ILC term from the last iteration

% Create new control input by adding perturbation to 
uNew = timesignal(timeseries(uNext,tscPath.pathVar.Data));



%% Plot some things
path        = lemOfGerono(linspace(0,1),basisParams);
tscTime.posVec.plot3('LineWidth',1,'Color','b','LineStyle','-','DisplayName','Flight Path')
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