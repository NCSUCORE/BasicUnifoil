clc

indx = 3;

[As,Bs] = pathLinearize(tscPath,basisParams,refGain1,refGain2,ca,cb,tauRef,baseMass+addedMass);
format compact
format long

phi     = tscPath.azimuth.Data(:);
theta   = tscPath.elevation.Data(:);
v       = tscPath.speed.Data(:);
psi     = tscPath.twistAngle.Data(:);
omega   = tscPath.twistRate.Data(:);
s       = tscPath.pathVar.Data(:);
u       = tscPath.twistSP.Data(:);
W       = basisParams(1);
H       = basisParams(2);
r       = basisParams(5);
a       = refGain1;
b       = refGain2;
c       = ca;
d       = cb;
tau     = tauRef;
M       = baseMass+addedMass;

varNames = {'phi','theta','v','psi','omega','s','u'};
mathmatExpr = '{';
for ii = 1:numel(varNames)
    strTmp = strcat(varNames{ii},sprintf('(%d)',indx));
    mathmatExpr = strcat(mathmatExpr,varNames{ii},'->',sprintf('%.9f',eval(strTmp)));
    mathmatExpr = strcat(mathmatExpr,',');
end

varNames = {'W','H','r','a','b','c','d','tau','M'};
for ii = 1:numel(varNames)
    strTmp = varNames{ii};
    mathmatExpr = strcat(mathmatExpr,varNames{ii},'->',sprintf('%.9f',eval(strTmp)));
    mathmatExpr = strcat(mathmatExpr,',');
end
mathmatExpr = mathmatExpr(1:end-1);
mathmatExpr = strcat(mathmatExpr,'};');
fprintf('\nSubstitution expression for Mathematica:\n')
fprintf(mathmatExpr)
fprintf('\n\nMatlab Answer')



format compact
format short
for ii = 1:5
    fprintf('\nColumn %d\n',ii)
    As.Data(1,ii,indx)
    As.Data(2,ii,indx)
    As.Data(3,ii,indx)
    As.Data(4,ii,indx)
    As.Data(5,ii,indx)
end

Bs.Data(:,:,indx)
