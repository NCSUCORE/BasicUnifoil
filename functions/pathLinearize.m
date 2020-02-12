function [Ats,Bts] = pathLinearize(tsc,basisParams,refGain1,refGain2,ca,cb,taur,mass)
% Rename a bunch of variables, otherwise this is intractable
phi     = tsc.azimuth.Data(:);
theta   = tsc.elevation.Data(:);
v       = tsc.speed.Data(:);
psi     = tsc.twistAngle.Data(:);
omega   = tsc.twistRate.Data(:);
s       = tsc.pathVar.Data(:);
psiSP   = tsc.twistSP.Data(:);
W       = basisParams(1);
H       = basisParams(2);
r       = basisParams(5);
ba      = refGain1;
bb      = refGain2;
M       = mass;

numSteps = numel(tsc.azimuth.Time);

% Intermediate variables, z's
za = W.*cos(2*pi*s).*cos(psi).*sec(theta);
zb = H.*cos(4*pi*s).*sin(phi);
zc = W^2.*cos(4*pi*s)+4.*H.*cos(8*pi*s) - 8*H*sin(4*pi*s).*theta+2*W*sin(2*pi*s).*phi;
zd = 8.*H.*sin(4.*pi.*s).*(-za+2.*zb)-za.*tan(theta).*(W.*cos(4.*pi.*s)+4.*H^2.*cos(8.*pi.*s)-8.*H.*sin(4.*pi.*s).*theta+2.*W.*sin(2.*pi.*s).*phi);
ze = W^2.*cos(2.*pi.*s).*cos(4.*pi.*s).*cos(phi).*cos(psi).*sec(theta)+...
    4.*H^2.*cos(2.*pi.*s).*cos(8.*pi.*s).*cos(phi).*cos(psi).*sec(theta)+...
    W.*cos(psi).*sec(theta).*sin(4.*pi.*s).*sin(phi)-...
    4.*H.*cos(4.*pi.*s).*sin(2.*pi.*s).*sin(phi).^2-...
    8.*H.*cos(2.*pi.*s).*cos(phi).*cos(psi).*sec(theta).*sin(4.*pi.*s).*theta+...
    W.*cos(phi).*cos(psi).*sec(theta).*sin(4.*pi.*s).*phi;
zf = 2.*H.*W^2.*cos(4.*pi.*s).^2.*cos(phi)+...
    8.*H^3.*cos(4.*pi.*s).*cos(8.*pi.*s).*cos(phi)+...
    W^2.*cos(psi).*sec(theta).*sin(4.*pi.*s)-...
    4.*W.*sin(2.*pi.*s).*zb-...
    8.*H^2.*cos(phi).*sin(8.*pi.*s).*theta+...
    4.*H.*W.*cos(4.*pi.*s).*cos(phi).*sin(2.*pi.*s).*phi;
zh = -4.*H.*sin(8.*pi.*s).*sin(phi)+...
    4.*sin(4.*pi.*s).*za+...
    W^2.*cos(4.*pi.*s).^2.*sin(phi).*tan(theta)+...
    4.*H^2.*cos(4.*pi.*s).*cos(8.*pi.*s).*sin(phi).*tan(theta)-...
    4.*H.*sin(8.*pi.*s).*sin(phi).*tan(theta).*theta+...
    2.*W.*cos(4.*pi.*s).*sin(2.*pi.*s).*sin(phi).*tan(theta).*phi;
zi = 1./((za-2.*zb).^2);


A = zeros(5,5,numSteps);
B = zeros(5,1,numSteps);

% Column 1
num = pi.*cos(psi).*sec(theta).*zf.*zi;
den = 1;
A(1,1,:) = num./den;

num = pi.*W.*ze.*zi;
den = 1;
A(2,1,:) = num./den;

num = exp(-(phi/cb)).*pi.*r.*cos(psi).^2.*ca.*(2.*W.*sin(2.*pi.*s).*(za-2.*zb)+2.*H.*cos(4.*pi.*s).*cos(phi).*zc-((za-2.*zb).*zc)/cb).*zi;
den = M.*v;
A(3,1,:) = num./den;

num = 2.*pi.*r.*(W.*sin(2.*pi.*s).*(za-2.*zb)+H.*cos(4.*pi.*s).*cos(phi).*zc).*zi.*omega;
den = v;
A(4,1,:) = num./den;

num = 2.*pi.*r.*(W.*sin(2.*pi.*s).*(za-2.*zb)+H.*cos(4.*pi.*s).*cos(phi).*zc).*zi.*(-ba.*psi-bb.*omega+psiSP/taur^2);
den = v;
A(5,1,:) = num./den;

% Column 2
num = -2*H*pi*cos(psi).*sec(theta).*zh.*zi;
den = 1;
A(1,2,:) = num./den;

num = pi*sin(phi).*zd.*zi;
den = 1;
A(2,2,:) = num./den;

num = exp(-phi./cb)*pi*r.*cos(psi).^2.*ca.*zd.*zi;
den = M.*v;
A(3,2,:) = num./den;

num = pi.*r.*zd.*zi.*omega;
den = v;
A(4,2,:) = num./den;

num = pi.*r.*zd.*zi.*(-ba.*psi-bb.*omega+psiSP./taur^2);
den = v;
A(5,2,:) = num./den;

% Column 3
% num = 0;
% den = 0;
% A(1,3,:) = num./den;
% 
% num = 0;
% den = 0;
% A(2,3,:) = num./den;

num = -exp(-phi./cb)*pi*r.*cos(psi).^2.*ca.*zc;
den = M.*(za-2.*zb).*v.^2;
A(3,3,:) = num./den;

num = -pi.*r.*zc.*omega;
den = (za-2.*zb).*v.^2;
A(4,3,:) = num./den;

num = -pi.*r.*zc.*(-ba.*psi-bb.*omega+psiSP./taur^2);
den = (za-2.*zb).*v.^2;
A(5,3,:) = num./den;

% Column 4
num = 2.*pi.*sec(theta).*sin(psi).*zb.*(W^2.*cos(4*pi*4)+4.*H^2.*cos(8*pi*s)-8.*H.*sin(4*pi*s).*theta+2.*W.*sin(2*pi*s).*phi);
den = (W.*cos(2*pi*s).*cos(psi).*sec(theta)-2.*H.*cos(4*pi*s).*sin(phi)).^2;
A(1,4,:) = num./den;

num = pi.*W.*cos(2*pi*s).*sec(theta).*sin(phi).*sin(psi).*zc.*zi;
den = 1;
A(2,4,:) = num./den;

num = exp(-phi/cb).*pi.*r.*ca.*(-W.*cos(2*pi*s).*cos(psi).^2.*sec(theta).*sin(psi)+2.*sin(psi).*zb).*zc.*zi;
den = M.*v;
A(3,4,:) = num./den;

num = pi.*r.*W.*cos(2*pi*s).*sec(theta).*sin(psi).*zc.*zi.*omega;
den = v;
A(4,4,:) = num./den;

num = -pi.*r.*zc.*zi.*(ba.*taur.^2.*(za-2.*zb+W.*cos(2.*pi.*s).*sec(theta).*sin(psi).*psi)+W.*cos(2.*pi.*s).*sec(theta).*sin(psi).*(bb.*taur.^2.*omega-psiSP));
den = v.*taur^2;
A(5,4,:) = num./den;

% Column 5
% num = 0;
% den = 0;
% A(1,5,:) = num./den;
% 
% num = 0;
% den = 0;
% A(2,5,:) = num./den;
% 
% num = 0;
% den = 0;
% A(3,5,:) = num./den;

num = pi.*r.*zc;
den = (za-2.*zb).*v;
A(4,5,:) = num./den;

num = -pi.*r.*bb.*zc;
den = (za-2.*zb).*v;
A(5,5,:) = num./den;

num = pi*r*(W^2*cos(4*pi*s) + 4*H^2*cos(8*pi*s) - 8*H*sin(4*pi*s).*theta + 2*W*sin(2*pi*s).*phi);
den = (W*cos(2.*pi.*s).*cos(psi).*sec(theta) - 2.*H.*cos(4.*pi.*s).*sin(phi)).*v.*taur^2;
B(5,1,:) = num./den;

Ats = timesignal(timeseries(A,tsc.stateVec.Time));
Bts = timesignal(timeseries(B,tsc.stateVec.Time));
end