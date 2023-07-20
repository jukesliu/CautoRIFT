#!/usr/bin/env python
# coding: utf-8

# In[1]:


def rio_write(out_path, nparray, ref_raster, grid_spacing):
    # Function to write a numpy array to a Geotiff using rasterio
    # The geotiff will have the same bounds and crs as the reference raster
    # An evenly-spaced grid will be created with the grid spacing entered
    #
    # INPUTS:
    #   out_path: the path with name of the output gtiff file
    #   nparray: the nparray to write to gtiff
    #   ref_raster: the reference raster (we borrow its crs and bounding coordinates)
    #   grid_spacing: the spatial resolution of the output raster
    
    import rasterio as rio
    import numpy as np
    nparray = np.array(nparray) # make sure it's an np array
    
    with rio.open(out_path,'w',
                      driver='GTiff',
                      height=nparray.shape[0], # new shape
                      width=nparray.shape[1], # new shape
                      dtype=nparray.dtype, # data type
                      count=1,
                      crs=ref_raster.crs, # the EPSG from the original DEM
                      transform=rio.Affine(grid_spacing, 0.0, ref_raster.bounds.left, # modified transform
                                           0.0, -grid_spacing, ref_raster.bounds.top)) as dst:
            dst.write(nparray, 1)


# In[2]:


def download_orbits(SAFEzipfilepath, config_path, out_dir):
    # DOWNLOADS PRECISE ORBIT FILES FROM ASF
    # Requires installation of wget
    # INPUTS:
    #  SAFEzipfilepath: path to the SAFE file 
    #  configpath: 
    #  our_dir: path to the orbit directory where orbit files will be saved
    import datetime
    import urllib.request
    import re
    import os
    import subprocess
    
    orb_type = 'aux_poeorb'
    zipname = SAFEzipfilepath.split('/')[-1] 
    time = zipname.split('_')[5]
    S1 = zipname.split('_')[0][-3:]
    scene_center_time = datetime.datetime.strptime(time,"%Y%m%dT%H%M%S")
    validity_start_time = scene_center_time-datetime.timedelta(days=1)
    validity_end_time =  scene_center_time+datetime.timedelta(days=1)
    
    # ASF URL
    url = "https://s1qc.asf.alaska.edu/%s/?validity_start=%s&validity_start=%s&validity_start=%s&sentinel1__mission=%s" % (orb_type, validity_start_time.strftime("%Y"),validity_start_time.strftime("%Y-%m"), validity_start_time.strftime("%Y-%m-%d"), S1)   
    content = ( urllib.request.urlopen(url).read()) # read results
    ii = re.findall('''href=["'](.[^"']+)["']''', content.decode('utf-8'))
    
    for i in ii :
        if '.EOF' in i:
            if (validity_start_time.strftime("%Y%m%d") in i) and (validity_end_time.strftime("%Y%m%d") in i) and (S1 in i):
                
                # if it doesn't already exist
                if not os.path.isfile(out_dir+i):
                    wget_cmd = 'export WGETRC="'+config_path+'"; '
                    wget_cmd += 'wget -c -P '+out_dir+' '
                    wget_cmd += "https://s1qc.asf.alaska.edu/aux_poeorb/"+i
    #                 print(wget_cmd)
                    subprocess.run(wget_cmd, shell=True,check=True)  
                    print(i+' downloaded.')
                else:
                    print(i+' already exists in orbit folder.')


# In[3]:


# create .netrc file with Earthdata credentials
def create_netrc(netrc_name):
    from netrc import netrc
    from subprocess import Popen
    from getpass import getpass
    import os
    
    homeDir = os.path.expanduser("~")
    if os.path.exists(homeDir+'/'+netrc_name):
        print(netrc_name, ' with Earthdata credentials already exists.')

    urs = 'urs.earthdata.nasa.gov'    # Earthdata URL endpoint for authentication
    prompts = ['Enter NASA Earthdata Login Username: ',
               'Enter NASA Earthdata Login Password: ']
    # Determine if netrc file exists, and if so, if it includes NASA Earthdata Login Credentials
    try:
        netrcDir = os.path.expanduser(f"~/{netrc_name}")
        netrc(netrcDir).authenticators(urs)[0]
    # Below, create a netrc file and prompt user for NASA Earthdata Login Username and Password
    except FileNotFoundError:
        Popen('touch {0}{2} | echo machine {1} >> {0}{2}'.format(homeDir + os.sep, urs, netrc_name), shell=True)
        Popen('echo login {} >> {}{}'.format(getpass(prompt=prompts[0]), homeDir + os.sep, netrc_name), shell=True)
        Popen('echo \'password {} \'>> {}{}'.format(getpass(prompt=prompts[1]), homeDir + os.sep, netrc_name), shell=True)
        # Set restrictive permissions
        Popen('chmod 0600 {0}{1}'.format(homeDir + os.sep, netrc_name), shell=True)


# In[ ]:

def generate_geogrid_inputs(CHIPSIZE_M, dempath, demname, refvpath, vx_fname, vy_fname, sr_scaling):
    
    # GRAB DEM INFO
    refdem = rio.open(dempath+demname) # open DEM using rasterio
    elev = refdem.read(1) # read in the first and only band (elevations)

    # grab the x and y grid values from the DEM:
    dem_x = np.linspace(refdem.bounds.left, refdem.bounds.right, num=np.shape(elev)[1])
    dem_y = np.linspace(refdem.bounds.top, refdem.bounds.bottom, num=np.shape(elev)[0])

    # grab the resampled x and y grid values from the DEM
    new_x = np.arange(refdem.bounds.left, refdem.bounds.right, CHIPSIZE_M)
    new_y = np.arange(refdem.bounds.top, refdem.bounds.bottom, -CHIPSIZE_M)
    
    # RESAMPLE THE DEM
    dem_outfile = 'IfSAR_'+str(CHIPSIZE_M)+'m_DSM_clipped.tif' # generate new filename
    if not os.path.exists(dempath+dem_outfile): # if the resampled DEM does not exist already
        # Create thew new x and y grid values using DEM bounds and the chipsize
        dem_resamp = np.zeros((len(new_y), len(new_x))) # create an empty resampled DEM grid
        print(dem_resamp.shape)

        # Resample to your new DEM bounds
        f = interp2d(dem_x, dem_y, elev) # create DEM interpolation object
        dem_resamp = f(new_x,new_y) # resample the NIR data to the DSM coordinates
        dem_resamp = np.flipud(dem_resamp) # flip up down
        print("Resampled to new dimensions:",dem_resamp.shape)

        # Display the two DEMs as a visual check
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5))
        im1 = ax1.imshow(elev, cmap='Greys_r', vmin=0)
        ax1.set_title('Original DEM: '+str(refdem.transform[0])+' m') # original spatial resolution
        fig.colorbar(im1, ax=ax1,label='Elevation [m]')

        im2 = ax2.imshow(dem_resamp, cmap='Greys_r', vmin=0)
        ax2.set_title('Resampled DEM: '+str(CHIPSIZE_M)+' m') # new spatial resolution
        fig.colorbar(im2, ax=ax2,label='Elevation [m]')
        plt.show()

        # Save the resampled DEM to georeferenced tif file
        print("Save resampled DEM to", dempath+dem_outfile)
        rio_write(dempath+dem_outfile, dem_resamp, refdem, CHIPSIZE_M)
    else:
        # load the existing resampled DEM
        dem_r = rio.open(dempath+dem_outfile) # open DEM using rasterio
        dem_resamp = dem_r.read(1) # read in the first and only band (elevations)
        print(dem_outfile, ' already exists.')
    
    # CREATE DHDX, DHDY
    dhdx_outfile = 'IfSAR_'+str(CHIPSIZE_M)+'m_DSM_clipped_dhdx.tif' # generate new filename
    dhdy_outfile = 'IfSAR_'+str(CHIPSIZE_M)+'m_DSM_clipped_dhdy.tif' # generate new filename
    if not os.path.exists(dempath+dhdx_outfile) or not os.path.exists(dempath+dhdy_outfile): # if either is missing
        # Produce dhdx and dhdy maps from resampled DEM
        dhdx = np.gradient(dem_resamp, axis=1)/CHIPSIZE_M
        dhdy = np.gradient(dem_resamp, axis=0)/CHIPSIZE_M

        # Filter out borders with high gradient values
        grad_thresh = 5
        dhdx[abs(dhdx) > grad_thresh] = 0; dhdy[abs(dhdy) > grad_thresh] = 0

        # absolute value of the max gradient values expected:
        dhmax = 1

        # Display the two DEMs as a visual check
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5))
        im1 = ax1.imshow(dhdx, cmap='Greys_r', vmin=-dhmax, vmax=dhmax)
        ax1.set_title('dhdx') # surface slope x
        fig.colorbar(im1, ax=ax1)

        im2 = ax2.imshow(dhdy, cmap='Greys_r', vmin=-dhmax, vmax=dhmax)
        ax2.set_title('dhdy') # surface slope y
        fig.colorbar(im2, ax=ax2)
        plt.show()

        # Save the gradient maps to tif files
        print("Save surface slope maps to", dempath)
        rio_write(dempath+dhdx_outfile, dhdx, refdem, CHIPSIZE_M) # dhdx
        rio_write(dempath+dhdy_outfile, dhdy, refdem, CHIPSIZE_M)
    else:
        print(dhdy_outfile, 'and', dhdx_outfile, 'already exist.')

    # VX, VY, SRX, SRY
    vx_outfile = 'vx_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    vy_outfile = 'vy_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    srx_outfile = 'srx_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    sry_outfile = 'sry_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    if not os.path.exists(refvpath+vx_outfile) or not os.path.exists(refvpath+vy_outfile): # if either vx, vy missing
        # open the files with rasterio
        vx_reader = rio.open(refvpath+vx_fname); vx0 = vx_reader.read(1)
        vy_reader = rio.open(refvpath+vy_fname); vy0 = vy_reader.read(1)
        vx_x = np.linspace(vx_reader.bounds.left, vx_reader.bounds.right, num=np.shape(vx0)[1])
        vx_y = np.linspace(vx_reader.bounds.top, vx_reader.bounds.bottom, num=np.shape(vx0)[0])
        vy_x = np.linspace(vy_reader.bounds.left, vy_reader.bounds.right, num=np.shape(vy0)[1])
        vy_y = np.linspace(vy_reader.bounds.top, vy_reader.bounds.bottom, num=np.shape(vy0)[0])

        # Resample to the DEM grid
        fx = interp2d(vx_x, vx_y, vx0)
        fy = interp2d(vy_x, vy_y, vy0)
        vx_resamp = np.flipud(fx(new_x,new_y)) 
        vy_resamp = np.flipud(fy(new_x,new_y)) # flip up down

        # Display the two velocity files as a visual check
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5))
        im1 = ax1.imshow(vx_resamp, cmap='Greys_r'); ax1.set_title('vx'); fig.colorbar(im1, ax=ax1)
        im2 = ax2.imshow(vy_resamp, cmap='Greys_r'); ax2.set_title('vy'); fig.colorbar(im2, ax=ax2)
        plt.show()

        # CALCULATE SEARCH RANGE LIMITS MULTIPLY VX AND VY BY SOME NUMBER
        srx_resamp = vx_resamp*sr_scaling; sry_resamp = vy_resamp*sr_scaling

        # Display the two search range files as a visual check
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5))
        im1 = ax1.imshow(srx_resamp, cmap='Greys_r'); ax1.set_title('srx'); fig.colorbar(im1, ax=ax1)
        im2 = ax2.imshow(sry_resamp, cmap='Greys_r'); ax2.set_title('sry'); fig.colorbar(im2, ax=ax2)
        plt.show()

        # save the reference velocity and search range maps
        rio_write(refvpath+vx_outfile, vx_resamp, refdem, CHIPSIZE_M) # vx
        rio_write(refvpath+vy_outfile, vy_resamp, refdem, CHIPSIZE_M) # vy
        rio_write(refvpath+srx_outfile, srx_resamp, refdem, CHIPSIZE_M) # srx
        rio_write(refvpath+sry_outfile, sry_resamp, refdem, CHIPSIZE_M) # sry
    else:
        print(vx_outfile, ',', vy_outfile, ',', srx_outfile, ',', sry_outfile, 'already exist.')  
    
    # MASKS
    ssm_outfile = 'ssm_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    tg_outfile = 'TG_mask_'+str(CHIPSIZE_M)+'m.tif' # generate new filename
    tg_mt_outfile = 'TG_mask_MT_'+str(CHIPSIZE_M)+'m.tif'
    tg_nt_outfile = 'TG_mask_NT_'+str(CHIPSIZE_M)+'m.tif'
    tg_st_outfile = 'TG_mask_ST_'+str(CHIPSIZE_M)+'m.tif'

    if True: # overwrite all
        # read it in, process (resample, mask, etc.) and resave
        ssmreader = rio.open(refvpath+'ice_mask_200mbuffer.tif')
        ssm = ssmreader.read(1)
        ssm[ssm > 0] = 1; #ssm[ssm < 0.0] = 0; # make binary
        ssm = ssm < 1 # find all stable areas (where.tif = 0)

        # do the same for Turner Glacier mask
        tgreader = rio.open(refvpath+'TG_mask.tif'); tg_mask = tgreader.read(1)
        tg_mask[tg_mask > 0] = 1; tg_mask = tg_mask > 0
        
        # and all the regional masks
        tg_mt_reader = rio.open(refvpath+'TG_mask_MT.tif'); tg_mt = tg_mt_reader.read(1)
        tg_mt[tg_mt > 0] = 1; tg_mt = tg_mt > 0
        tg_nt_reader = rio.open(refvpath+'TG_mask_NT.tif'); tg_nt = tg_nt_reader.read(1)
        tg_nt[tg_nt > 0] = 1; tg_nt = tg_nt > 0
        tg_st_reader = rio.open(refvpath+'TG_mask_ST.tif'); tg_st = tg_st_reader.read(1)
        tg_st[tg_st > 0] = 1; tg_st = tg_st > 0
        
#         fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(10,5)); ax1.imshow(tg_mt); ax1.set_title('main truk')
#         ax2.imshow(tg_nt); ax2.set_title('N tributary'); ax3.imshow(tg_st); ax3.set_title('S tributary')
#         plt.show()

        # grab x and y-values
        ssm_x = np.linspace(ssmreader.bounds.left, ssmreader.bounds.right, num=np.shape(ssm)[1])
        ssm_y = np.linspace(ssmreader.bounds.top, ssmreader.bounds.bottom, num=np.shape(ssm)[0])
        tg_x = np.linspace(tgreader.bounds.left, tgreader.bounds.right, num=np.shape(tg_mask)[1])
        tg_y = np.linspace(tgreader.bounds.top, tgreader.bounds.bottom, num=np.shape(tg_mask)[0])
        tg_mtx = np.linspace(tg_mt_reader.bounds.left, tg_mt_reader.bounds.right, num=np.shape(tg_mt)[1])
        tg_mty = np.linspace(tg_mt_reader.bounds.top, tg_mt_reader.bounds.bottom, num=np.shape(tg_mt)[0])

        # Resample to the DEM grid
        f_ssm = interp2d(ssm_x, ssm_y, ssm)
        f_tg = interp2d(tg_x, tg_y, tg_mask)
        f_tg_mt = interp2d(tg_mtx, tg_mty, tg_mt)
        f_tg_nt = interp2d(tg_mtx, tg_mty, tg_nt)
        f_tg_st = interp2d(tg_mtx, tg_mty, tg_st)
        
        ssm_resamp = np.flipud(f_ssm(new_x,new_y))
        tg_resamp = np.flipud(f_tg(new_x, new_y))
        tg_mt_resamp = np.flipud(f_tg_mt(new_x, new_y))
        tg_nt_resamp = np.flipud(f_tg_nt(new_x, new_y))
        tg_st_resamp = np.flipud(f_tg_st(new_x, new_y))
        
        # plot
        fig, ax = plt.subplots(1,1)
        ssm_im = ax.imshow(ssm_resamp,cmap='gray',vmin=0)
        ax.set_title('Stable Surface Mask')
        fig.colorbar(ssm_im, ax=ax)
        plt.show()

        # export
        rio_write(refvpath+ssm_outfile, ssm_resamp, refdem, CHIPSIZE_M)
        rio_write(refvpath+tg_outfile, tg_resamp, refdem, CHIPSIZE_M)
        rio_write(refvpath+tg_mt_outfile, tg_mt_resamp, refdem, CHIPSIZE_M)
        rio_write(refvpath+tg_nt_outfile, tg_nt_resamp, refdem, CHIPSIZE_M)
        rio_write(refvpath+tg_st_outfile, tg_st_resamp, refdem, CHIPSIZE_M)
    else:
        print(ssm_outfile,'and',tg_outfile,'already exist.')
    
    return dem_outfile, dhdx_outfile, dhdy_outfile, vx_outfile, vy_outfile, srx_outfile, sry_outfile, ssm_outfile, tg_outfile
        
def run_geogrid_inhouse(out_path, img_type, indir_m, indir_s, MINCHIPSIZE, NO_DATA_VAL, dem, # required inputs
                        dhdx, dhdy, vx, vy, srx, sry, csminx, csminy, csmaxx, csmaxy, ssm, # optional inputs
                       temp_dir): # SAR needed only
    
    CHIPSIZE_M = MINCHIPSIZE # set minimum chip size equal
    ############ Clear all old geogrid files ##########################
    for file in os.listdir(out_path):
        if file.startswith('window') and file.endswith('.tif'):
            print('removed', file)
            os.remove(out_path+file)
    print('Old files cleared.'); print()

    dem_info = gdal.Info(dem, format='json') # grab info from DEM
    print('Obtained DEM info.'); print()

    ############ Run geogrid optical or SAR ##########################
    if img_type == 'OPT': # Optical images
        print('Processing optical images with geogrid.'); print()
        obj = GeogridOptical() # initialize geogrid object

        ############ Coregister the optical data (from coregisterLoadMetadataOptical) #############
        x1a, y1a, xsize1, ysize1, x2a, y2a, xsize2, ysize2, trans = obj.coregister(indir_m, indir_s,0)

        # grab dates from file names
        im1_name = indir_m.split('/')[-1]; im2_name = indir_s.split('/')[-1]
        if 'LC' in im1_name and 'LC' in im2_name:
            ds1 = im1_name.split('_')[3]
            ds2 = im1_name.split('_')[4]
        elif 'S2' in im1_name and 'S2' in im2_name:
            ds1 = im1_name[9:17]
            ds2 = im2_name[9:17]
        elif 'PS' in im1_name and 'PS' in im2_name:
            ds1 = im1_name[3:11]
            ds2 = im2_name[3:11]
        else:
            raise Exception('Optical data NOT supported yet!') 
        print('Optical images coregistered.'); print()

        ########### Load geogrid inputs and run (from runGeogridOptical) ################

        # grab info from above
        obj.startingX = trans[0]; obj.startingY = trans[3]
        obj.XSize = trans[1]; obj.YSize = trans[5]
        d0 = datetime.date(int(ds1[0:4]),int(ds1[4:6]),int(ds1[6:8]))
        d1 = datetime.date(int(ds2[0:4]),int(ds2[4:6]),int(ds2[6:8]))
        date_dt_base = d1 - d0
        obj.repeatTime = date_dt_base.total_seconds()
        obj.numberOfLines = ysize1; obj.numberOfSamples = xsize1
        obj.gridSpacingX = dem_info['geoTransform'][1] # output grid spacing is the same as the DEM

        # customize no data value and minimimum chip size
        obj.nodata_out = NO_DATA_VAL
        obj.chipSizeX0 = MINCHIPSIZE

        # set raster paths and names
        obj.dat1name = indir_m # first image
        obj.demname = dem # DEM
        obj.dhdxname = dhdx; obj.dhdyname = dhdy # surface slope
        obj.vxname = vx; obj.vyname = vy # reference velocity
        obj.srxname = srx; obj.sryname = sry # search range limits
        obj.csminxname = csminx; obj.csminyname = csminy # min chip size
        obj.csmaxxname = csmaxx; obj.csmaxyname = csmaxy # max chip size
        obj.ssmname = ssm # stable surface mask
        obj.winlocname = "window_location.tif"
        obj.winoffname = "window_offset.tif"
        obj.winsrname = "window_search_range.tif"
        obj.wincsminname = "window_chip_size_min.tif"
        obj.wincsmaxname = "window_chip_size_max.tif"
        obj.winssmname = "window_stable_surface_mask.tif"
        obj.winro2vxname = "window_rdr_off2vel_x_vec.tif"
        obj.winro2vyname = "window_rdr_off2vel_y_vec.tif"

        obj.runGeogrid() # RUN GEOGRID
        print('Optical geogrid finished.'); print()

    elif img_type == 'SAR': # SAR images
        print('Processing SAR images with geogrid.'); print();
        ############ Load SAR metadata from coreg_files ##################################
        # Store sensing start info for 2nd SAR image (in temp_dir+secondary/)
        frames = []
        for swath in range(1,4):
            inxml = os.path.join(temp_dir+'secondary/', 'IW{0}.xml'.format(swath))
            if os.path.exists(inxml):
                pm = PM(); pm.configure(); ifg = pm.loadProduct(inxml) # load XML file
                frames.append(ifg)
        info1_sensingStart = min([x.sensingStart for x in frames]) # store info1_sensingStart

        # Load other info from 1st SAR image (in temp_dir+reference/)
        del frames; frames = [] 
        for swath in range(1,4):
            inxml = os.path.join(temp_dir+'reference/', 'IW{0}.xml'.format(swath))
            if os.path.exists(inxml):
                pm = PM(); pm.configure(); ifg = pm.loadProduct(inxml) # load XML file        
                frames.append(ifg)
        print('SAR metadata loaded.'); print()

        ############ Get merged orbit getMergedOrbit() ################################## 
        # Create merged orbit
        orb = Orbit(); orb.configure()
        burst = frames[0].bursts[0]
        # Add first burst orbit to begin with
        for sv in burst.orbit:
            orb.addStateVector(sv)
        for pp in frames:
            # Add all state vectors
            for bb in pp.bursts:
                for sv in bb.orbit:
                    if (sv.time< orb.minTime) or (sv.time > orb.maxTime):
                        orb.addStateVector(sv)
        print('Merged orbit created.'); print()

        ############ Load geogrid inputs and run ###################################
        obj = Geogrid()
        obj.configure()

        obj.orbit = orb # grab merged orbit
        obj.startingRange = min([x.startingRange for x in frames])
        obj.rangePixelSize = frames[0].bursts[0].rangePixelSize
        obj.sensingStart = min([x.sensingStart for x in frames])
        obj.prf = 1.0 / frames[0].bursts[0].azimuthTimeInterval
        obj.lookSide = -1
        obj.repeatTime = (info1_sensingStart - obj.sensingStart).total_seconds() # INFO1
        obj.numberOfLines = int(np.round((max([x.sensingStop for x in frames])-obj.sensingStart).total_seconds()*obj.prf))+1
        obj.numberOfSamples = int(np.round((max([x.farRange for x in frames])-obj.startingRange)/obj.rangePixelSize))+1
        obj.gridSpacingX = dem_info['geoTransform'][1] # output grid spacing is the same as the DEM

        # custom no data value and chip size
        obj.nodata_out = NO_DATA_VAL
        obj.chipSizeX0 = CHIPSIZE_M

        # set raster paths and names
        obj.demname = dem # DEM
        obj.dhdxname = dhdx; obj.dhdyname = dhdy # surface slope
        obj.vxname = vx; obj.vyname = vy # reference velocity
        obj.srxname = srx; obj.sryname = sry # search range limmits
        obj.csminxname = csminx; obj.csminyname = csminy # min chip size
        obj.csmaxxname = csmaxx; obj.csmaxyname = csmaxy # max chip size
        obj.ssmname = ssm # stable surface mask
        obj.winlocname = "window_location.tif"
        obj.winoffname = "window_offset.tif"
        obj.winsrname = "window_search_range.tif"
        obj.wincsminname = "window_chip_size_min.tif"
        obj.wincsmaxname = "window_chip_size_max.tif"
        obj.winssmname = "window_stable_surface_mask.tif"
        obj.winro2vxname = "window_rdr_off2vel_x_vec.tif"
        obj.winro2vyname = "window_rdr_off2vel_y_vec.tif"

        obj.getIncidenceAngle() # SAR specific
        obj.geogrid() # run geogrid
        print('SAR geogrid finished.'); print();

    else: # not OPT or SAR
        print('Image type flag not recognized :', img_type)


    ############ Move files produced to the out_path directory ##############
    for file in os.listdir(os.getcwd()):
        if file.startswith('window') and file.endswith('.tif'):
            shutil.move(os.getcwd()+'/'+file, out_path+file)
    print('Geogrid output files moved')
    

# AutoRIFT  
def run_autoRIFT_inhouse(out_path, img_type, mpflag, xGrid, yGrid, indir_m, indir_s, # required parameters
                         FILTER, WALLISFILTERWIDTH, SPARSE_SEARCH_SAMPLE_RATIO, OVERSAMPLE_RATIO, MINCHIPSIZE,
                         Dx0, Dy0, CSMINx0, SRx0, SRy0, CSMAXx0, CSMAXy0, SSM, # optional parameters
                         noDataMask, nodataval, geogrid_run_info):
    import logging
    from imageMath import IML
    
    
    CHIPSIZE_M = MINCHIPSIZE # set minimum chip size equal
    
    # requires grid location from geogrid
    origSize = xGrid.shape # grab original size from xGrid
    
    if img_type == 'OPT': ############# OPTICAL SETTINGS ############################# 
        print('Processing optical images with autoRIFT.'); print()
        optflag = 1 # turn on optical flag
        # Coregister and read in the two images (from loadProductOptical())
        obj = GeogridOptical()
        x1a, y1a, xsize1, ysize1, x2a, y2a, xsize2, ysize2, trans = obj.coregister(indir_m, indir_s,0)

        # read dates from filenames
        if 'LC' in indir_m and 'LC' in indir_s:
            ds1 = indir_m.split('/')[-1].split('_')[3]; ds2 = indir_s.split('/')[-1].split('_')[3]
            sat = 'LS'
        elif 'S2' in indir_m and 'S2' in indir_s:
            ds1 = indir_m.split('/')[-1][9:17]; ds2 = indir_s.split('/')[-1][9:17]
            sat = 'S2'
        elif 'PS' in im1_name and 'PS' in im2_name:
            ds1 = indir_m.split('/')[-1].split('_')[1]
            ds2 = indir_s.split('/')[-1].split('_')[1]
            sat = 'PS'
        else:
            raise Exception('Optical data NOT supported yet!')

        # read in the images
        DS1 = gdal.Open(indir_m); DS2 = gdal.Open(indir_s)
        I1 = DS1.ReadAsArray(xoff=x1a, yoff=y1a, xsize=xsize1, ysize=ysize1)
        I1 = I1.astype(np.float32)
        I2 = DS2.ReadAsArray(xoff=x2a, yoff=y2a, xsize=xsize2, ysize=ysize2)
        I2 = I2.astype(np.float32)
        DS1=None; DS2=None # clear DS1 and DS2

    elif img_type == 'SAR': ############# SAR SETTINGS #############################  
        print('Processing SAR images with autoRIFT.'); print()
        optflag = 0 # turn off opt flag
        # Read in the two SAR images (from loadProduct())
        img1 = IML.mmapFromISCE(indir_m, logging); I1 = img1.bands[0]
        img2 = IML.mmapFromISCE(indir_s, logging); I2 = img2.bands[0]
        I1 = np.abs(I1); I2 = np.abs(I2) # SAR amplitude only
        ds1 = indir_m.split('/')[-1].split('.')[0]
        ds2 = indir_s.split('/')[-1].split('.')[0]
        sat = 'S1'
    else:
        print("Image type not recognized. Use either 'OPT' or 'SAR'.")
        
    ############# Initialize autoRIFT object (from runAutorift()) ##################
    obj = autoRIFT_ISCE()
    obj.configure()
    
    obj.MultiThread = mpflag # multiprocessing
    obj.I1 = I1; obj.I2 = I2 # assign the images
    obj.xGrid = xGrid; obj.yGrid = yGrid # assign the grid 

    # GENERATE NO DATA MASK
    # where offset searching will be skipped based on 
    # 1) imported nodata mask and/or 2) zero values in the image
    for ii in range(obj.xGrid.shape[0]):
        for jj in range(obj.xGrid.shape[1]):
            if (obj.yGrid[ii,jj] != nodata)&(obj.xGrid[ii,jj] != nodata):
                if (I1[obj.yGrid[ii,jj]-1,obj.xGrid[ii,jj]-1]==0)|(I2[obj.yGrid[ii,jj]-1,obj.xGrid[ii,jj]-1]==0):
                    noDataMask[ii,jj] = True
                    
    # SEARCH RANGE
    if SRx0 is None:
        # default is a zero array
#        ###########     uncomment to customize SearchLimit based on velocity distribution 
        if Dx0 is not None:
            obj.SearchLimitX = np.int32(4+(25-4)/(np.max(np.abs(Dx0[np.logical_not(noDataMask)]))-np.min(np.abs(Dx0[np.logical_not(noDataMask)])))*(np.abs(Dx0)-np.min(np.abs(Dx0[np.logical_not(noDataMask)]))))
        else:
            obj.SearchLimitX = 15
        obj.SearchLimitY = 15
#        ###########
        obj.SearchLimitX = obj.SearchLimitX * np.logical_not(noDataMask)
        obj.SearchLimitY = obj.SearchLimitY * np.logical_not(noDataMask)
    else:
        obj.SearchLimitX = SRx0
        obj.SearchLimitY = SRy0
       ############ add buffer to search range
        obj.SearchLimitX[obj.SearchLimitX!=0] = obj.SearchLimitX[obj.SearchLimitX!=0] + 2
        obj.SearchLimitY[obj.SearchLimitY!=0] = obj.SearchLimitY[obj.SearchLimitY!=0] + 2
    
    # CHIP SIZE
    if CSMINx0 is not None:
        obj.ChipSizeMaxX = CSMAXx0
        obj.ChipSizeMinX = CSMINx0
        
        gridspacingx = MINCHIPSIZE # use the grid spacing from above
        if img_type == 'OPT':
            pixsizex = trans[1] # grab from coregister function
            chipsizex0 = MINCHIPSIZE
        else:
            pixsizex = 2.3 # 2.3 m single look IW pixel spacing in rnage direction
            chipsizex0 = np.unique(CSMINx0)[1]
            
        obj.ChipSize0X = int(np.ceil(chipsizex0/pixsizex/4)*4)
        obj.GridSpacingX = int(obj.ChipSize0X*gridspacingx/chipsizex0)

        RATIO_Y2X = CSMINy0/CSMINx0
        obj.ScaleChipSizeY = np.median(RATIO_Y2X[(CSMINx0!=nodata)&(CSMINy0!=nodata)])
#         obj.ScaleChipSizeY = 1 # USE SCALE OF 1 for square pixels
    else:
        if ((optflag == 1)&(xGrid is not None)):
            obj.ChipSizeMaxX = 32 # pixels
            obj.ChipSizeMinX = 16 # pixels
            obj.ChipSize0X = 16 # pixels
    
    # DOWNSTREAM SEARCH OFFSET
    if Dx0 is not None:
        obj.Dx0 = Dx0
        obj.Dy0 = Dy0
    else:
        obj.Dx0 = obj.Dx0 * np.logical_not(noDataMask)
        obj.Dy0 = obj.Dy0 * np.logical_not(noDataMask)

    # REPLACE NO DATA VALUES WITH 0
    obj.xGrid[noDataMask] = 0
    obj.yGrid[noDataMask] = 0
    obj.Dx0[noDataMask] = 0
    obj.Dy0[noDataMask] = 0
    if SRx0 is not None:
        obj.SearchLimitX[noDataMask] = 0
        obj.SearchLimitY[noDataMask] = 0
    if CSMINx0 is not None:
        obj.ChipSizeMaxX[noDataMask] = 0
        obj.ChipSizeMinX[noDataMask] = 0
    
    # convert azimuth offset to vertical offset as used in autoRIFT convention for SAR images
    if optflag == 0:
        obj.Dy0 = -1 * obj.Dy0
        
    ############## AutoRIFT Pre-processing (from runAutorift()) ############################
    t1 = time.time()
    print("Pre-process Start!!!")
    
    # FILTERING:
    if FILTER == 'WAL': 
        obj.preprocess_filt_wal() # WALLIS FILTER
#         obj.zeroMask = 1 # removes edges
        obj.WallisFilterWidth = WALLISFILTERWIDTH # optional, default supposedly 21
    elif FILTER == 'HPS':
        obj.preprocess_filt_hps() # HIGH PASS FILTER
    elif FILTER == 'SOB':
        obj.preprocess_filt_sob() # SOBEL FILTER
    elif FILTER == 'LAP':
        obj.preprocess_filt_lap()
    elif FILTER == 'DB':
        obj.preprocess_db() # LOGARITHMIC OPERATOR (NO FILTER), FOR TOPOGRAPHY
    else:
        print(FILTER, 'not recognized. Using default high pass filter instead.')
        obj.preprocess_filt_hps() # HIGH PASS FILTER
        
    print("Pre-process Done!!!")
    print(time.time()-t1)
    
    # CONVERT TO UNIFORM DATA TYPE
    t1 = time.time()
#    obj.DataType = 0
    obj.uniform_data_type()
    print("Uniform Data Type Done!!!")
    print(time.time()-t1)
    
    # OTHER :
    obj.sparseSearchSampleRate = 1
#    obj.colfiltChunkSize = 4

    obj.OverSampleRatio = 64
    if CSMINx0 is not None:
        obj.OverSampleRatio = {obj.ChipSize0X:16,obj.ChipSize0X*2:32,obj.ChipSize0X*4:64,obj.ChipSize0X*8:64}
    
    # SEE ORIGINAL CODE TO EXPORT PREPROCESSED IMAGES
    
    ####################### Run AutoRIFT (from runAutorift())  ############################
    t1 = time.time()
    print("AutoRIFT Start!!!")
    obj.runAutorift()
    print("AutoRIFT Done!!!")
    print(time.time()-t1)
    
    kernel = np.ones((3,3),np.uint8)
    noDataMask = cv2.dilate(noDataMask.astype(np.uint8),kernel,iterations = 1)
    noDataMask = noDataMask.astype(np.bool)

    # AT THIS POINT, THESE VARIABLES WILL BE CREATED:
    # obj.Dx, obj.Dy, obj.InterpMask, obj.ChipSizeX, obj.GridSpacingX, 
    # obj.ScaleChipSizeY, obj.SearchLimitX, obj.SearchLimitY, obj.origSize, noDataMask
    
    # PLOT RESULTS
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15,5))
    im1 = ax1.imshow(obj.Dx); ax1.set_title('Dx'); fig.colorbar(im1, ax=ax1)
    im2 = ax2.imshow(obj.Dy); ax2.set_title('Dy'); fig.colorbar(im2, ax=ax2)
    im3 = ax3.imshow(np.sqrt((obj.Dx**2) + (obj.Dy**2))); ax3.set_title('D_total'); fig.colorbar(im3,ax=ax3)
    plt.suptitle(ds1+' to '+ds2)
    plt.show()

    ####################### Write outputs (from runAutorift())  ############################
    t1 = time.time()
    print("Write Outputs Start!!!")
          
    # Write text file with parameters
    f =  open(out_path+'parameters_'+ds1+'_'+ds2+'_'+str(CHIPSIZE_M)+'m_'+sat+'.txt', 'w')
    f.write('Geogrid/AutoRIFT parameters for offset_'+ds1+'_'+ds2+'_'+str(CHIPSIZE_M)+'m_'+sat+'.tif:')
    f.write('NO_DATA_VAL: '+str(NO_DATA_VAL))
    f.write('Min chip size: '+str(MINCHIPSIZE))
    f.write('DEM: '+dem)
    f.write('dhdx: '+dhdx); f.write('dhdy: '+dhdy)
    f.write('vx: '+vx); f.write('vy: '+vy)
    f.write('srx: '+srx); f.write('sry: '+sry)
    f.write('csminx: '+csminx); f.write('csminy: '+csminy)
    f.write('csmaxx: '+csmaxx); f.write('csmaxy: '+csmaxy)
    f.write('stable surface mask: '+ssm)
    f.write('FILTER: '+FILTER)
    f.write('WALLISFILTERWIDTH: '+str(WALLISFILTERWIDTH))
    f.write('Spare search sample rate: '+str(SPARSE_SEARCH_SAMPLE_RATE))
    f.write('Oversample ratio: '+str(OVERSAMPLE_RATIO))
    if offset2vx is not None and offset2vy is not None:
        f.write('Velocity.TIF file created.')
    else:
        f.write('Velocity.TIF not created.')
    f.close() # close the parameter text file
          
    # open the window_location.tif file to gdalinfo
    ds = gdal.Open(gp+'window_location.tif')
    tran = ds.GetGeoTransform()
    proj = ds.GetProjection()
    srs = ds.GetSpatialRef()
    
    # initialize arrays
    DX = np.zeros(origSize,dtype=np.float32) * np.nan; DY = np.zeros(origSize,dtype=np.float32) * np.nan
    INTERPMASK = np.zeros(origSize,dtype=np.float32); CHIPSIZEX = np.zeros(origSize,dtype=np.float32)
    SEARCHLIMITX = np.zeros(origSize,dtype=np.float32); SEARCHLIMITY = np.zeros(origSize,dtype=np.float32)
    
    # fill in arays
    Dx = obj.Dx; Dy = obj.Dy; InterpMask = obj.InterpMask; ChipSizeX = obj.ChipSizeX
    SearchLimitX = obj.SearchLimitX; SearchLimitY = obj.SearchLimitY
    DX[0:Dx.shape[0],0:Dx.shape[1]] = Dx;  DY[0:Dy.shape[0],0:Dy.shape[1]] = Dy
    INTERPMASK[0:InterpMask.shape[0],0:InterpMask.shape[1]] = InterpMask
    CHIPSIZEX[0:ChipSizeX.shape[0],0:ChipSizeX.shape[1]] = ChipSizeX
    SEARCHLIMITX[0:SearchLimitX.shape[0],0:SearchLimitX.shape[1]] = SearchLimitX
    SEARCHLIMITY[0:SearchLimitY.shape[0],0:SearchLimitY.shape[1]] = SearchLimitY
    
    # mask out no data
    DX[noDataMask] = np.nan; DY[noDataMask] = np.nan
    INTERPMASK[noDataMask] = 0; CHIPSIZEX[noDataMask] = 0
    SEARCHLIMITX[noDataMask] = 0; SEARCHLIMITY[noDataMask] = 0
    if SSM is not None:
        SSM[noDataMask] = False
    DX[SEARCHLIMITX == 0] = np.nan; DY[SEARCHLIMITX == 0] = np.nan
    INTERPMASK[SEARCHLIMITX == 0] = 0; CHIPSIZEX[SEARCHLIMITX == 0] = 0
    if SSM is not None:
        SSM[SEARCHLIMITX == 0] = False

    # SAVE TO OFFSET.MAT FILE
    sio.savemat('offset_'+ds1+'_'+ds2+'_'+str(CHIPSIZE_M)+'m_'+sat+'.mat', # offset mat filename
                {'Dx':DX,'Dy':DY,'InterpMask':INTERPMASK,'ChipSizeX':CHIPSIZEX})
    print('Offset.mat written.')
    
    # CREATE THE GEOTIFFS
    driver = gdal.GetDriverByName('GTiff')
    
    # OFFSET.TIF
    outRaster = driver.Create("offset_"+ds1+'_'+ds2+'_'+str(CHIPSIZE_M)+'m_'+sat+".tif", # offset filename
                              int(xGrid.shape[1]), int(xGrid.shape[0]), 5, gdal.GDT_Float32)
    outRaster.SetGeoTransform(tran); outRaster.SetProjection(proj) # projections
    outband = outRaster.GetRasterBand(1); outband.WriteArray(DX) # DX
    outband.FlushCache()
    outband = outRaster.GetRasterBand(2); outband.WriteArray(DY) # DY
    outband.FlushCache()
    outband = outRaster.GetRasterBand(3); outband.WriteArray(np.sqrt((DX**2) + (DY**2))) # DY
    outband.FlushCache()
    outband = outRaster.GetRasterBand(4); outband.WriteArray(INTERPMASK) # INTERPMASK
    outband.FlushCache()
    outband = outRaster.GetRasterBand(5); outband.WriteArray(CHIPSIZEX) # CHIPSIZE
    outband.FlushCache()
    del outRaster
    print('Offset.tif written.')
    
    # VELOCITY.TIF
    if offset2vx is not None and offset2vy is not None:
        ds = gdal.Open(offset2vx) #### VX
        band = ds.GetRasterBand(1); offset2vx_1 = band.ReadAsArray()
        band = ds.GetRasterBand(2); offset2vx_2 = band.ReadAsArray()
        if ds.RasterCount > 2:
                band = ds.GetRasterBand(3)
                offset2vr = band.ReadAsArray()
        else:
                offset2vr = None
        band=None; ds=None
        offset2vx_1[offset2vx_1 == nodata] = np.nan
        offset2vx_2[offset2vx_2 == nodata] = np.nan

        ds = gdal.Open(offset2vy) #### VY
        band = ds.GetRasterBand(1); offset2vy_1 = band.ReadAsArray()
        band = ds.GetRasterBand(2); offset2vy_2 = band.ReadAsArray()
        if ds.RasterCount > 2:
                band = ds.GetRasterBand(3)
                offset2va = band.ReadAsArray()
        else:
                offset2va = None
        band=None; ds=None
        offset2vy_1[offset2vy_1 == nodata] = np.nan; offset2vy_2[offset2vy_2 == nodata] = np.nan
        
        if offset2va is not None:
            offset2va[offset2va == nodata] = np.nan

        VX = offset2vx_1 * DX + offset2vx_2 * DY
        VY = offset2vy_1 * DX + offset2vy_2 * DY
        VX = VX.astype(np.float32); VY = VY.astype(np.float32)

        outRaster = driver.Create("velocity_"+ds1+'_'+ds2+'_'+str(CHIPSIZE_M)+'m_'+sat+".tif", # velocity filename
                                  int(xGrid.shape[1]), int(xGrid.shape[0]), 3, gdal.GDT_Float32)
        outRaster.SetGeoTransform(tran); outRaster.SetProjection(proj)
        outband = outRaster.GetRasterBand(1); outband.WriteArray(VX) # VX
        outband.FlushCache()
        outband = outRaster.GetRasterBand(2); outband.WriteArray(VY) # VY
        outband.FlushCache()
        outband = outRaster.GetRasterBand(3); outband.WriteArray(np.sqrt((VX**2) + (VY**2))) # V
        outband.FlushCache()
        del outRaster
        print('Velocity.tif written.')
    
    print("Write Outputs Done!!!")
    print(time.time()-t1)
    
    # Move files produced to the out_path directory
    for file in os.listdir(os.getcwd()):
        if 'offset' in file or ('velocity' in file and '.tif' in file):
            shutil.move(os.getcwd()+'/'+file, out_path+file)

