function gndPos = lemOfGerono(pathVar,geomParams)
%LEMOFGERONO Mercator projection of the lemniscate of Gerono onto a sphere
%   INPUTS
%   pathVar = Normalized path variabl (0-1) describing position along the
%   path
%   geomParams(1) = A0, total azimuth sweep angle in degrees
%   geomParams(2) = Z0, total elevation sweep angle in degrees
%   geomParams(3) = A1, mean course azimuth angle in degrees
%   geomParams(4) = Z1, mean course elevation angle in degrees
%   geomParams(5) = R, radius of sphere
%   OUTPUTS
%   gndPos = Nx3 matrix, where N is the number of elements of pathVar
%   Variable names here were chosen to match the notebook
%   lemOfGeronoTanVec.nb

A0 = geomParams(1)*pi/180;
Z0 = geomParams(2)*pi/180;
A1 = geomParams(3)*pi/180;
Z1 = geomParams(4)*pi/180;
R  = geomParams(5);

% Calculate path position in rad, phi from the normalized path variable, s.
phi = (pathVar(:)*2+3/2)*pi;

% Calculate azimuth and zenith in degrees
a =         (A0/2)*cos(  phi(:)) + A1;
z = (pi/2)-((Z0/2)*sin(2*phi(:)) + Z1);

% Convert sphereical to cartesian
% http://mathworld.wolfram.com/SphericalCoordinates.html
gndPos = nan(numel(phi),3);

gndPos(:,1) = R.*cos(a).*sin(z);
gndPos(:,2) = R.*sin(a).*sin(z);
gndPos(:,3) = R.*cos(z);

end
