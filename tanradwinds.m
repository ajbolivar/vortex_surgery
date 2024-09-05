function [utan,urad]=tanradwinds(u,v)
% tanradwinds.m - Written by Jake Carstens 9/7/2022
% Takes u and v wind components from shear-relative composite scripts, and converts
% to tangential and radial winds. This is called at each relevant vertical level
% individually from a 10 degree-wide TC-centered box.
% Outputs tangential and radial winds on the same grid as the u/v inputs.

  [nx,ny]=size(u); xc=floor(nx/2)+1; yc=floor(ny/2)+1; % Establish grid size/TC center
  [xx,yy] = meshgrid(1:nx,1:ny); theta = atan2(xx-xc,yy-yc);

% Calculating radial and tangential winds on cartesian grid
  urad = u.*cos(theta)+v.*sin(theta); % radial wind
  utan = -u.*sin(theta)+v.*cos(theta); % tangential wind
end
