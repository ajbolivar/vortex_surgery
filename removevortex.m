% Program RemoveVortex_ERA5_TEST.m
% Originally written by Jake Carstens 3/7/2023, using vortex flow removal code by Nick Barron and Chris Yu
% Functionized by Ana Bolivar 4/22/2024 to be used with a python wrapper.

% Purpose: Remove vortex flow from tropical cyclones in a reanalysis,
%          by removing the non-divergent and irrotational winds in a TC-centered annulus
%          whose radius is prescribed by a radius of vanishing axisymmetric tangential winds.
%          This process yields a vortex-free wind field representative of the synoptic
%          environment, which after possible adjustments (i.e. smoothing functions, tweaks to
%          the box size), will allow 6-hourly climate data to be used in a beta advection model
%          for TC tracks, rather than less representative monthly means.

% STEP 1: Load in center fixes from trajectories_ERA5.txt derived from TempestExtremes
% STEP 2: For each center fix, load in corresponding winds at 850/200 hPa and calculate tangential component
% STEP 3: Remove vortex flow using helmholtz_decompose.m, courtesy of Nick Barron
% STEP 4: Smooth wind field where vortex was removed (not yet implemented).
% STEP 5: Output .nc file with modified wind fields

function [] = removevortex(vr,xgrid,ygrid,radius,ufile,vfile,u850name,v850name,u200name,v200name,lonname,latname,i,j,t,outfile)
% --Inputs:
%   xgrid: longitude grid spacing
%   ygrid: latitude grid spacing
%   radius: radius of smoothing around storm center (km)
%   ufile: filepath of u wind data (if u and v are stored together, set ufile and vfile to the same path)
%   vfile: filepath of v wind data (if u and v are stored together, set ufile and vfile to the same path)
%   u850name: variable name for 850 hPa u winds
%   v850name: variable name for 850 hPa v winds
%   u200name: variable name for 200 hPa u winds
%   v200name: variable name for 200 hPa v winds
%   lonname: variable name for longitude
%   latname: variable name for latitude
%   i: storm center latitude index
%   j: storm center longitude index
%   t: time index
%
% First version: 15 Jan 2012
% Updated: 22 Apr 2024
%--------------------------------------------------------------------------
xbox=7.5./xgrid; ybox=7.5./ygrid;
U=ncread(ufile,u850name); 
V=ncread(vfile,v850name); 
U850=U(:,:,t); 
V850=V(:,:,t); 
U850new=U850; 
V850new=V850;
U=ncread(ufile,u200name); V=ncread(vfile,v200name); U200=U(:,:,t); V200=V(:,:,t); U200new=U200; V200new=V200;
times=ncread(ufile,'time'); time=times(t); 
time_units=ncreadatt(ufile,'time','units'); time_cal=ncreadatt(ufile,'time','calendar')
lat_start=ncread(ufile,latname); lon_start=ncread(ufile,lonname);
xtotal=length(lon_start)
nstorms=length(i);
clear U V

%if vortex removal is true, proceed with calculations
if vr == 1
    for n=1:nstorms
      tic
      xint=[]; yint=[]; lat_at_time=[]; lon_at_time=[];
      if (i(n) < xbox+1)
        xstart=xtotal-(xbox-i(n));
      else
        xstart=i(n)-xbox;
      end
      if (i(n) > xtotal-xbox)
        xend=xbox-(xtotal-i(n));
      else
        xend=i(n)+xbox;
      end
      yint=j(n)-ybox:j(n)+ybox;
      if (xstart > xtotal-(2*xbox))
        xint=[xstart:xtotal 1:xend];
      else
        xint=xstart:xend;
      end
    
      for x=1:length(xint)
        lat_at_time(x,:)=lat_start(yint);
      end
      for x=1:length(yint)
        lon_at_time(:,x)=lon_start(xint);
      end
    
      %calculate the storm centers.
      centerfixlat = lat_at_time(floor(length(xint)/2)+1,floor(length(yint)/2)+1);
      centerfixlon = lon_at_time(floor(length(xint)/2)+1,floor(length(yint)/2)+1);
      %put the storm centers into an array that matches the size of the wind data.
      centerfixlonmesh = ones(size(lon_at_time))*centerfixlon;
      centerfixlatmesh = ones(size(lat_at_time))*centerfixlat;
    
      %calculate the distance from each point from the center of the storm.
      [xd, yd, ~,~] = latlon_to_disaz(centerfixlatmesh, centerfixlonmesh, lat_at_time, lon_at_time);
    
      dis = sqrt(xd.^2+yd.^2);
    
      %mask to find all distances less than the defined smoothing radius
      mask = dis/1000<radius;
    
      [dx, dy, ~,~] = latlon_to_disaz(lat_at_time(1,1), lon_at_time(1,1), lat_at_time(2,2), lon_at_time(2,2));
    
      %calculate the helmholtz decomposition of the winds at each level
      [u850_psi,v850_psi,u850_chi,v850_chi,psi850,chi850]=helmholtz_decompose(U850(xint,yint),V850(xint,yint),dx,dy);
      [u200_psi,v200_psi,u200_chi,v200_chi,psi200,chi200]=helmholtz_decompose(U200(xint,yint),V200(xint,yint),dx,dy);
    
      for a = 1:length(xint)
        for b = 1:length(yint)
            if mask(a, b) == 1
                % Calculate distance from center of TC for each point
                distance = lldistkm([lat_start(yint(b)) lon_start(xint(a))], [centerfixlat centerfixlon]);
                % Calculate the sigmoid coefficient
                k = 0.005;
                d0 = radius / 2; % Midpoint
                coeff = 1 / (1 + exp(-k * (d0 - distance)));
    
                U850new(xint(a), yint(b)) = coeff * (U850(xint(a), yint(b)) - u850_psi(a, b) - u850_chi(a, b)) + (1 - coeff) * U850(xint(a), yint(b));
                V850new(xint(a), yint(b)) = coeff * (V850(xint(a), yint(b)) - v850_psi(a, b) - v850_chi(a, b)) + (1 - coeff) * V850(xint(a), yint(b));
                U200new(xint(a), yint(b)) = coeff * (U200(xint(a), yint(b)) - u200_psi(a, b) - u200_chi(a, b)) + (1 - coeff) * U200(xint(a), yint(b));
                V200new(xint(a), yint(b)) = coeff * (V200(xint(a), yint(b)) - v200_psi(a, b) - v200_chi(a, b)) + (1 - coeff) * V200(xint(a), yint(b));
    
                clear distance, coeff
            end
        end
      end      
    end
    
    nccreate(outfile, 'lon', 'Dimensions', {'lon', length(lon_start)});
    nccreate(outfile, 'lat', 'Dimensions', {'lat', length(lat_start)});
    nccreate(outfile, 'time', 'Dimensions', {'time', 0});
    nccreate(outfile, u850name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, v850name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, u200name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, v200name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    
    ncwriteatt(outfile, 'time', 'units', time_units);
    ncwriteatt(outfile, 'time', 'calendar', time_cal);
    
    ncwrite(outfile, u850name, U850new);
    ncwrite(outfile, v850name, V850new);
    ncwrite(outfile, u200name, U200new);
    ncwrite(outfile, v200name, V200new);
    ncwrite(outfile, 'lon', lon_start);
    ncwrite(outfile, 'lat', lat_start);
    ncwrite(outfile, 'time', time)
%if vortex removal is false, output the original file
else
    nccreate(outfile, 'lon', 'Dimensions', {'lon', length(lon_start)});
    nccreate(outfile, 'lat', 'Dimensions', {'lat', length(lat_start)});
    nccreate(outfile, 'time', 'Dimensions', {'time', 0});
    nccreate(outfile, u850name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, v850name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, u200name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    nccreate(outfile, v200name, 'Dimensions', {'lon', length(lon_start), 'lat', length(lat_start), 'time', 0});
    
    ncwriteatt(outfile, 'time', 'units', time_units);
    ncwriteatt(outfile, 'time', 'calendar', time_cal);
    
    ncwrite(outfile, u850name, U850);
    ncwrite(outfile, v850name, V850);
    ncwrite(outfile, u200name, U200);
    ncwrite(outfile, v200name, V200);
    ncwrite(outfile, 'lon', lon_start);
    ncwrite(outfile, 'lat', lat_start);
    ncwrite(outfile, 'time', time)
end
