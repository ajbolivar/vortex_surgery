vr = 0;
xgrid = 0.25;
ygrid = 0.25;
radius = 700;
ufile = "/glade/derecho/scratch/abolivar/h1files/1989/ERA5.h1.19890426.nc";
vfile = "/glade/derecho/scratch/abolivar/h1files/1989/ERA5.h1.19890426.nc";
u850name = "U850";
v850name = "V850";
u200name = "U200";
v200name = "V200";
lonname = "lon";
latname = "lat";
i = -999;
j = -999;
t = 1;
outfile = "/glade/derecho/scratch/abolivar/h1files/vortex_removed/1989/ERA5.h1.1989042600.nc.tmp"
removevortex(vr,xgrid,ygrid,radius,ufile,vfile,u850name,v850name,u200name,v200name,lonname,latname,i,j,t,outfile)