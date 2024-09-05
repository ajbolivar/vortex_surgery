import numpy as np
import pandas as pd
import xarray as xr
import subprocess
import os
import sys

# Function that creates .m file with variables to be used in removevortex.m
def vortex_removal(mfn, vr, xgrid, ygrid, radius, ufn, vfn, u850name, v850name, 
                   u200name, v200name, lonname, latname, i, j, t, outfile):
    m_script = open(mfn, 'w')
    script  = f'vr = {vr};\n'
    script  += f'xgrid = {xgrid};\n'
    script  += f'ygrid = {ygrid};\n'
    script  += f'radius = {radius};\n'
    script  += f'ufile = "{ufn}";\n'
    script  += f'vfile = "{vfn}";\n'
    script  += f'u850name = "{u850name}";\n'
    script  += f'v850name = "{v850name}";\n'
    script  += f'u200name = "{u200name}";\n'
    script  += f'v200name = "{v200name}";\n'
    script  += f'lonname = "{lonname}";\n'
    script  += f'latname = "{latname}";\n'
    script  += f'i = {i};\n'
    script  += f'j = {j};\n'
    script  += f't = {t};\n'
    script  += f'outfile = "{outfile}"\n'
    script  += 'removevortex(vr,xgrid,ygrid,radius,ufile,vfile,u850name,v850name,u200name,v200name,lonname,latname,i,j,t,outfile)'
    m_script.write(script)
    m_script.close()

    sout = subprocess.run(['bash', 'run_vortex_removal.sh'],
                           stdout=subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines=True)
    return sout

# Arguments for matlab script
xgrid = 0.25
ygrid = 0.25
radius = 700
u850name = 'U850'
v850name = 'V850'
u200name = 'U200'
v200name = 'V200'
lonname = 'lon'
latname = 'lat'
mfn = 'vortex_removal_args.m'
traj_file = '/glade/work/abolivar/Pyclogenesis_data/track_data/ERA5/ERA5_trajectories_1979-2014.txt'
ufile = '/glade/derecho/scratch/abolivar/h1files/{}/ERA5.h1.{}{:02d}{:02d}.nc'
vfile = '/glade/derecho/scratch/abolivar/h1files/{}/ERA5.h1.{}{:02d}{:02d}.nc'
outfile = '/glade/derecho/scratch/abolivar/h1files/vortex_removed/{}/ERA5.h1.{}{:02d}{:02d}{:02d}.nc'

year_to_process = int(sys.argv[1])  # Get the year from the command line argument

print("Building trajectories file...")
# Build trajectories file
td = open(traj_file)
lines = td.readlines()
# Split lines by tab and remove new line characters
lines_split = [l.replace('\n','').split('\t')for l in lines]
# Find all start line indexes and append last line too
lines_start = [i for i in range(0,len(lines_split)) if lines_split[i][0] == 'start'] + [len(lines_split)] 
storm_length = list(np.diff(lines_start) - 1) # Minus one because start line would be included if not
trajectories = pd.read_table(traj_file, skiprows = lines_start[:-1], usecols = ([3,4,8,9,10,11]),
                            delimiter = '\t', header=None, 
                            names = ['lon','lat','yyyy','mm','dd','hh'])

# Create time array
time_array = (trajectories['yyyy']*1000000+
              trajectories['mm']*10000+
              trajectories['dd']*100+
              trajectories['hh']*1)

trajectories['time'] = pd.to_datetime(time_array, format='%Y%m%d%H')
# Make sure all lons and lats are on a 0.25 degree grid (for edge cases)
trajectories['lon']  = round(trajectories['lon'] * int(1 / xgrid)) / (1 / xgrid)
trajectories['lat']  = round(trajectories['lat'] * int(1 / ygrid)) / (1 / ygrid)
# Remove year, month, day, and hour columns
trajectories = trajectories.drop(['yyyy','mm','dd','hh'], axis=1)
trajectories = trajectories.drop_duplicates()

# Subset trajectories by year
trajectories = trajectories.loc[trajectories['time'].dt.year == year_to_process]
# Create array of times for the whole year
dates = pd.date_range(f'{year_to_process}-01-01', f'{year_to_process + 1}-01-01', freq = '6H', inclusive = 'left')
# dates = pd.date_range('1985-06-26', '1985-06-27', freq = '6H')
# Loop over all times, performing vortex removal on timesteps with vortices present
for date in dates:
    print(f'Current date/time: {date}')
    year = pd.to_datetime(date).year
    month = pd.to_datetime(date).month
    day = pd.to_datetime(date).day
    hour = pd.to_datetime(date).hour

    ufn = ufile.format(year,year,month,day)
    vfn = vfile.format(year,year,month,day)
    ofn = outfile.format(year,year,month,day,hour)
    temp_ofn = ofn + ".tmp"
    
    # Skip files that already exist
    if os.path.exists(ofn):
        print(f'File {ofn} already exists, skipping...')
        continue

    os.makedirs(f'/glade/derecho/scratch/abolivar/h1files/vortex_removed/{year}', exist_ok=True)
        
    # Skip dates if they don't exist
    if os.path.exists(ufn) and os.path.exists(vfn):
        u = xr.open_dataset(ufn)
        lon  = u.lon
        lat  = u.lat
        time = u.time
        
        # Add 1 to indices for matlab
        t = np.argwhere(time.values==date).item() + 1

        if date in trajectories:
            vr = 1
            traj_at_time = trajectories.loc[trajectories['time'] == date]
            if len(traj_at_time) > 1:
                i = list(np.squeeze([np.argwhere(lon.values==traj_lon)+1 for traj_lon in traj_at_time.lon]))
                j = list(np.squeeze([np.argwhere(lat.values==traj_lat)+1 for traj_lat in traj_at_time.lat]))
            elif len(traj_at_time) == 1:
                i = [np.argwhere(lon.values==traj_at_time.lon.item()).item()+1]
                j = [np.argwhere(lat.values==traj_at_time.lat.item()).item()+1]

        else:
            vr = 0
            # No vortex removal will occur, so pass dummy values
            i = -999
            j = -999
        
        # Perform the operation and check if the file was saved correctly
        success = False
        attempts = 0
        while not success:
            if attempts >= 10: break
            attempts += 1
            sout = vortex_removal(mfn, vr, xgrid, ygrid, radius, ufn, vfn, u850name, v850name, u200name, v200name, lonname, latname, i, j, t, temp_ofn) 
            # print(sout.stderr)

            # Use ncdump to check if the temp file is valid
            result = subprocess.run(['ncdump', '-h', temp_ofn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:  # ncdump was successful
                success = True
                os.rename(temp_ofn, ofn)  # Move temp file to final destination
            else:
                print(f"File generation failed for {ofn}, retrying...")
                try:
                    os.remove(temp_ofn)  # Remove the corrupted file before retrying
                except:
                    continue

        if success: print(f"File generation succeeded for {ofn}!")
        elif not success: print(f"Too many attempts for {ofn}. Check source file for corruption.")
