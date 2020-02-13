function out = dsdtinv(tsc,basisParams)
phi     = tsc.azimuth.Data(:);
theta   = tsc.elevation.Data(:);
v       = tsc.speed.Data(:);
psi     = tsc.twistAngle.Data(:);
s       = tsc.pathVar.Data(:);
W       = basisParams(1);
H       = basisParams(2);
r       = basisParams(5);

num = (W*cos(2*pi*s).*cos(psi).*sec(theta)-2*H*cos(4*pi*s).*sin(phi)).*v;
den = pi*r*(W^2*cos(4*pi*s)+4*H^2*cos(8*pi*s)-8*H*sin(4*pi*s).*theta+2*W*sin(2*pi*s).*phi);

out = den./num;
end