function [Ats,Bts] = pathLinearize(tsc,basisParams,FxParams,tauRef,mass)
% The math behind this block is derived in the matlab notebook
% Reparameterization.nb.  This notebook automatically outputs .txt files
% containing matlab expressions for each element of A and B.

% Rename variables to match variable names used in mathematica
phi     = tsc.azimuth.Data(:);
theta   = tsc.elevation.Data(:);
v       = tsc.speed.Data(:);
psi     = tsc.twistAngle.Data(:);
omega   = tsc.twistRate.Data(:);
s       = tsc.pathVar.Data(:);
u       = tsc.twistSP.Data(:);
tau     = tauRef;
M       = mass;

a = FxParams(1);
b = FxParams(2);
c = basisParams(1);
d = basisParams(3);
e = basisParams(2);
f = basisParams(4);
r = basisParams(5);

shapeSwitch = basisParams(6);
if shapeSwitch == 1
    fldr = 'Ellipse';
else
    fldr = 'Fig8';
end

% Preallocate A and B matrices
A = zeros(5,5,numel(tsc.pathVar.Time));
B = zeros(5,1,numel(tsc.pathVar.Time));

% Find the files with the A and B matrices in them as calculated in
% Mathematica
basePath = fullfile(fileparts(which('unifoil.prj')),'Documentation');

% Evaluate all the expressions for the elements of A
for ii = 1:5
    for jj = 1:5
        fName = sprintf('A%d%d.txt',ii,jj);
        expr  = fileread(fullfile(basePath,'MatrixElementExpressions',fldr,fName));
        A(ii,jj,:) = eval(expr);
    end
end

% Evaluate all the expressions for the elements of B
for ii = 1:5
    fName = sprintf('B%d.txt',ii);
    expr  = fileread(fullfile(basePath,'MatrixElementExpressions',fldr,fName));
    B(ii,:) = eval(expr);
end

% Store everything into a timesignal for output
Ats = timesignal(timeseries(A,tsc.stateVec.Time));
Bts = timesignal(timeseries(B,tsc.stateVec.Time));
end
