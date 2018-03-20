import aplpy
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import pdb


#modify according the filters in the survey
def get_path_for_analysis_folder(gal,band,type_of_analysis,number_cpus):
    base_path = "/mnt/disk1/fb/gama_compacts/galaxies_all_selections_lee_18jul17/"

    #band
    if band == "g":
        path = base_path+"analysis_filter_g/"
    elif band == "r":
        path = base_path+"analysis_filter_r/"
    elif band == "i":
        path = base_path+"analysis_filter_i/"
    elif band == "z":
        path = base_path+"analysis_filter_z/"

    #type of analysis
    if type_of_analysis == "single":
        path = path + "single_sersic_fits/"
    elif type_of_analysis == "bulge_disk":
        path = path + "bulge_disk_decompositions/"

    for ii in range(1,number_cpus+1):
        total_path = path+"cpu"+str(ii)+"/"+str(gal)+"/"
        if os.path.exists(total_path):
            return(total_path)

#modify according to the stars in the survey
def get_stars_for_analysis_folder(gal,band):
    if band == "g":
        stars = np.array(["psf_for_"+str(gal)+"_g","star_in_185134_g","star_in_9204_g","star_in_23869_g"],dtype=str) #same order as they appear in autofit
    elif band == "r":
        stars = np.array(["psf_for_"+str(gal)+"_r","star_in_279707_r","star_in_138954_r","star_in_196099_r"],dtype=str) #same order as they appear in autofit
    elif band == "i":
        stars = np.array(["psf_for_"+str(gal)+"_i","star_in_23869_i","star_in_23880_i","star_in_388187_i"],dtype=str) #same order as they appear in autofit
    elif band == "z":
        stars = np.array(["star_in_371990_z","star_in_583311_z","star_in_740221_z"],dtype=str) #same order as they appear in autofit

    return(stars)

#modify according to the masks used in the survey
def get_masks_for_analysis_folder():
    masks = np.array(["normal_mask"],dtype=str) #same order as they appear in autofit
    return(masks)

def get_table_file_with_struct_param(band,type_of_analysis):
    path = "./analysis_filter_"+band+"/"
    if band != "z":
        if type_of_analysis == "single":
            file_to_read = "gama_kids_"+band+"_single_sersic_fits.cat"
        else:
            file_to_read = "gama_kids_"+band+"_bulge_disk_decompositions.cat"
    else:
        if type_of_analysis == "single":
            file_to_read = "viking_z_single_sersic_fits.cat"
        else:
            file_to_read = "viking_z_bulge_disk_decompositions.cat"

    return(path+file_to_read)

def get_filename(gal,band,type_of_analysis):
    stars = get_stars_for_analysis_folder(gal,band)
    masks = get_masks_for_analysis_folder()

    table_to_read = get_table_file_with_struct_param(band,type_of_analysis)

    if type_of_analysis == "single":
        tt = Table.read(table_to_read,
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index_of_best_fit"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            pass
        elif not np.isnan(index):
            #star_number==========================
            star_number = index #assuming no central_mask
            star_number = int(star_number)
            #=====================================
            #mask_number==========================
            if masks.size == 1:
                mask_number = 0
            else:
                mask_number = index % 2
            mask_number = int(mask_number)
            #=====================================
            #initial_conditions===================
            ini_cond_number = 0
            ini_cond_number = int(ini_cond_number)
            #=====================================
    elif type_of_analysis == "bulge_disk":
        tt = Table.read(table_to_read,
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center1","y_center1","mag1","mag_error1","mag1_exist","re_pix1",\
                          "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
                          "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","x_center2","y_center2","mag2",\
                          "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
                          "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
                          "pa_error2","pa2_exist","index_of_best_fit","flag_obj_detected"],
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index_of_best_fit"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            pass
        elif not np.isnan(index):
            #star_number==========================
            star_number = index / 2 #assuming no central_mask
            star_number = int(star_number)
            #=====================================
            #mask_number==========================
            if masks.size == 1:
                mask_number = 0
            else:
                pass #mask_number     = (index / 2) / 2
            mask_number = int(mask_number)
            #=====================================
            #initial_conditions===================
            ini_cond_number = index % 2
            ini_cond_number = int(ini_cond_number)
            #=====================================

    if np.isnan(index): #if the there is no best fit, return nothing
        return("")
    else:
        return( str(gal)+"_"+stars[star_number]+"_"+masks[mask_number]+"_"+str(ini_cond_number)+".fits" )

def create_sb_image(image,zp,pix_scale):
    sb = -2.5*np.log10(image)+zp+5.*np.log10(pix_scale)
    return(sb)

def select_conditions(band,camera,survey,special):

    zp_kids = 30.
    pix_scale_kids = 0.21
    fiveSigma_2arcsec_mag_limit_kids_g = 25.4 #from http://kids.strw.leidenuniv.nl/techspecs.php
    fiveSigma_2arcsec_mag_limit_kids_r = 25.2 #from http://kids.strw.leidenuniv.nl/techspecs.php
    fiveSigma_2arcsec_mag_limit_kids_i = 24.2 #from http://kids.strw.leidenuniv.nl/techspecs.php
    zp_viking_z = 30.
    pix_scale_viking_z = 0.339
    fiveSigma_mag_limit_viking_z = 23.1 #Venemans et al. (2013) Table 1

    if   band == "g":
        zp = zp_kids
        pix_scale = pix_scale_kids
        noise_level = fiveSigma_2arcsec_mag_limit_kids_g + 2.5*np.log10(5.)
    elif band == "r":
        zp = zp_kids
        pix_scale = pix_scale_kids
        noise_level = fiveSigma_2arcsec_mag_limit_kids_r + 2.5*np.log10(5.)
    elif band == "i":
        zp = zp_kids
        pix_scale = pix_scale_kids
        noise_level = fiveSigma_2arcsec_mag_limit_kids_i + 2.5*np.log10(5.)
    elif band == "z":
        zp = zp_viking_z
        pix_scale = pix_scale_viking_z
        noise_level = fiveSigma_mag_limit_viking_z + 2.5*np.log10(5.)
        noise_level = noise_level+2 #because the previous limit seems not to be enough

    return( (zp,pix_scale,noise_level) )

def get_struct_param_best_fit(gal,band,type_of_analysis):
    stars = get_stars_for_analysis_folder(gal,band)
    masks = get_masks_for_analysis_folder()

    table_to_read = get_table_file_with_struct_param(band,type_of_analysis)

    if type_of_analysis == "single":
        tt = Table.read(table_to_read,
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index_of_best_fit"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            return()
        else:
            index = int(index)
            re_pix = tt["re_pix"][filter_gal]; re_pix = np.float(re_pix[0])
            re_kpc = tt["re_kpc"][filter_gal]; re_kpc = np.float(re_kpc[0])
            nn     = tt["n"]     [filter_gal]; nn     = np.float(nn[0])
            ar     = tt["ar"]    [filter_gal]; ar     = np.float(ar[0])
            pa     = tt["pa"]    [filter_gal]; pa     = np.float(pa[0])
            return(re_pix,re_kpc,nn,ar,pa)
    elif type_of_analysis == "bulge_disk":
        tt = Table.read(table_to_read,
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center1","y_center1","mag1","mag_error1","mag1_exist","re_pix1",\
                          "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
                          "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","x_center2","y_center2","mag2",\
                          "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
                          "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
                          "pa_error2","pa2_exist","index_of_best_fit","flag_obj_detected"],
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index_of_best_fit"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            return()
        else:
            index = int(index)
            re1_pix = tt["re_pix1"][filter_gal]; re1_pix = np.float(re1_pix[0])
            re1_kpc = tt["re_kpc1"][filter_gal]; re1_kpc = np.float(re1_kpc[0])
            nn1     = tt["n1"]     [filter_gal]; nn1     = np.float(nn1[0])
            ar1     = tt["ar1"]    [filter_gal]; ar1     = np.float(ar1[0])
            pa1     = tt["pa1"]    [filter_gal]; pa1     = np.float(pa1[0])
            re2_pix = tt["re_pix2"][filter_gal]; re2_pix = np.float(re2_pix[0])
            re2_kpc = tt["re_kpc2"][filter_gal]; re2_kpc = np.float(re2_kpc[0])
            nn2     = tt["n2"]     [filter_gal]; nn2     = np.float(nn2[0])
            ar2     = tt["ar2"]    [filter_gal]; ar2     = np.float(ar2[0])
            pa2     = tt["pa2"]    [filter_gal]; pa2     = np.float(pa2[0])
            return(re1_pix,re1_kpc,nn1,ar1,pa1,re2_pix,re2_kpc,nn2,ar2,pa2)

def plotting_options():
    #-----axis-----
    subfig.axis_labels.hide()
    subfig.tick_labels.hide()
    subfig.ticks.hide()

def plotting_data_original(gal,band,filename,zz,logmass):
    subfig.add_label(0.50,0.95,"{0}"            .format(str(filename[:-5])),relative=True,color="black",weight="bold")
    subfig.add_label(0.70,0.85,"{0}"            .format(str(gal))          ,relative=True,color="black",weight="bold")
    subfig.add_label(0.70,0.80,"Filter: {0}"    .format(str(band))         ,relative=True,color="black",weight="bold")
    subfig.add_label(0.70,0.75,"z = {0:4.2f}"   .format(zz)                ,relative=True,color="black",weight="bold")
    subfig.add_label(0.70,0.70,"mass = {0:6.2e}".format(10.**logmass)      ,relative=True,color="black",weight="bold")
    
def plotting_data_single_model(re_kpc,re_pix,nn,ar):
    subfig.add_label(0.70,0.90,"re_kpc = {0:4.2f}"     .format(re_kpc)            ,relative=True,color='black',weight="bold")
    subfig.add_label(0.70,0.85,"re_pix = {0:4.2f}"     .format(re_pix)            ,relative=True,color='black',weight="bold")
    subfig.add_label(0.70,0.80,"re_circ_kpc = {0:4.2f}".format(re_kpc*np.sqrt(ar)),relative=True,color='black',weight="bold")
    subfig.add_label(0.70,0.75,"n      = {0:4.2f}"     .format(nn)                ,relative=True,color='black',weight="bold")

def plotting_data_double_model(re1_kpc,re1_pix,nn1,ar1,re2_kpc,re2_pix,nn2,ar2):
    subfig.add_label(0.80,0.90,r"r$_{e,bulge}$ = "+"{0:4.2f} kpc"     .format(re1_kpc)             ,relative=True,color='black')
    subfig.add_label(0.80,0.85,r"r$_{e,bulge}$ = "+"{0:4.2f} pix"     .format(re1_pix)             ,relative=True,color='black')
    subfig.add_label(0.80,0.80,r"r$_{e,circ,bulge}$ = "+"{0:4.2f} kpc".format(re1_kpc*np.sqrt(ar1)),relative=True,color='black')
    subfig.add_label(0.80,0.75,r"n$_{bulge}$ = "+"{0:4.2f}"           .format(nn1)                 ,relative=True,color='black')
    subfig.add_label(0.80,0.70,r"r$_{e,disk}$ = "+"{0:4.2f} kpc"      .format(re2_kpc)             ,relative=True,color='black')
    subfig.add_label(0.80,0.65,r"r$_{e,disk}$ = "+"{0:4.2f} pix"      .format(re2_pix)             ,relative=True,color='black')
    subfig.add_label(0.80,0.60,r"r$_{e,circ,disk}$ = "+"{0:4.2f} kpc" .format(re2_kpc*np.sqrt(ar2)),relative=True,color='black')
    subfig.add_label(0.80,0.55,r"n$_{disk}$ = "+"{0:4.2f}"            .format(nn2)                 ,relative=True,color='black')

def finding_coo_brightest_SB(image):
    #if the brightest pixel was the central, this function will only be return np.unravel_index(np.nanargmin( image ),image.shape)
    x_size = image[0].header["NAXIS1"]
    y_size = image[0].header["NAXIS2"]
    
    matrix = np.array(image[0].data)
    #looking in the inner part of the stamp
    offset_x = x_size/50
    offset_y = y_size/50
    submatrix = matrix[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
    indices_min_value = np.unravel_index(np.nanargmin(submatrix),submatrix.shape)
    filter_pos = matrix == submatrix[indices_min_value[0],indices_min_value[1]] #I assume that there is only a single pixel with the brightest value
    pos = np.where(filter_pos)
    return(pos[0][0]+1,pos[1][0]+1)
    #np.nanargmin->gets the minimum element of the image without taking into account the NaN values
    #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector

def printing_effective_radius(image,pix_scale,re_arcsec,ar,pa):
    offset_x_by_eye = 0. #(-1.*pix_scale)/3600.
    offset_y_by_eye = 0. #(1.*pix_scale)/3600.

    yy,xx=finding_coo_brightest_SB(image)
    #ellipses
    ra, dec = subfig.pixel2world(xx,yy)
    aa = re_arcsec/3600.
    bb = (re_arcsec*ar)/3600.
    pa = pa-90.
    #     pa=   catalog['pa_sex'][(ii*8)+(jj*4)]-90.
    subfig.show_ellipses(ra+offset_x_by_eye, dec+offset_y_by_eye, 2.*aa, 2.*bb, pa, edgecolor='#FFA500', linewidth=1)

def plotting_colorbar(box_position,max_value=None,min_value=None):
    subfig.add_colorbar()
    if max_value is not None:
        subfig.colorbar.show(box=box_position,box_orientation='horizontal',axis_label_text=r"Surf. brightness / mag arcsec$^{-2}$", ticks=np.round(np.linspace(max_value,min_value,6),decimals=1))
    else:
        subfig.colorbar.show(box=box_position,box_orientation='horizontal',axis_label_text=r"Surf. brightness ")   
    subfig.colorbar.set_font(weight="bold")
    subfig.colorbar.set_axis_label_font(weight="bold")

def plotting_scalebar():
    subfig.add_scalebar(1./3600.)
    subfig.scalebar.set_font(weight="bold")
    subfig.scalebar.set_color('black')
    subfig.scalebar.set_label('1 arcsec')
    
def str2bool(ele):
    return ele in ("True")


#CONSTANTS
catalog_path = "/mnt/disk1/fb/gama_compacts/galaxy_selection/"
catalog      = "galaxies_with_new_mass.cat"
bands = ["g","r","i","z"]
type_of_analysis = "single" #"bulge_disk" 
number_cpus = 3

plotting_pos= [[0.00,0.75,0.25,0.25],
               [0.25,0.75,0.25,0.25],
               [0.50,0.75,0.25,0.25],
               [0.75,0.75,0.25,0.25],
               [0.00,0.50,0.25,0.25],
               [0.25,0.50,0.25,0.25],
               [0.50,0.50,0.25,0.25],
               [0.75,0.50,0.25,0.25],
               [0.00,0.25,0.25,0.25],
               [0.25,0.25,0.25,0.25],
               [0.50,0.25,0.25,0.25],
               [0.75,0.25,0.25,0.25],
               [0.00,0.00,0.25,0.25],
               [0.25,0.00,0.25,0.25],
               [0.50,0.00,0.25,0.25],
               [0.75,0.00,0.25,0.25],
              ]


#reading the catalog
tt = Table.read(catalog_path+catalog, names=("CATAID","RA","DEC","Z_TONRY","logmstar","logmstar_Ferreras"), format='ascii.commented_header')
cataid  = tt["CATAID"]
ra      = tt["RA"]
dec     = tt["DEC"]
zz      = tt["Z_TONRY"]
logmass = tt["logmstar_Ferreras"]

#to see whether the fits were good or wrong
tt_good_or_bad = Table.read("final_sample_good_or_bad_fits.cat", \
                            names=("gal_id","final_good_fit_g","final_good_fit_r","final_good_fit_i","final_good_fit_z"), \
                            format='ascii.commented_header')

for cont_gal,gal in enumerate(cataid):

    fig = plt.figure( figsize = (16,16) )
    print(gal) #for debugging purposes

    #if the galaxy is also present in the final sample
    flag_final_sample = False
    filter_good_or_bad = tt_good_or_bad["gal_id"] == gal
    if filter_good_or_bad.any() == True:
        vector_good_fit = [str2bool(tt_good_or_bad["final_good_fit_g"][filter_good_or_bad][0]), str2bool(tt_good_or_bad["final_good_fit_r"][filter_good_or_bad][0]), str2bool(tt_good_or_bad["final_good_fit_i"][filter_good_or_bad][0]), str2bool(tt_good_or_bad["final_good_fit_z"][filter_good_or_bad][0])]
        filter_big_cat = cataid == gal
        flag_final_sample = True
        print(float(ra[filter_big_cat]),float(dec[filter_big_cat]),vector_good_fit)
    else:
        continue

    cont = 0
    for cont_band,band in enumerate(bands):

        path = get_path_for_analysis_folder(gal,band,type_of_analysis,number_cpus)

        stars = get_stars_for_analysis_folder(gal,band)

        filename = get_filename(gal,band,type_of_analysis)
        if filename == "": continue

        zp, pix_scale, noise_level = select_conditions(band,"","","")

        img = fits.open(path+filename)
        original_image  = img[1].data
        original_header = img[1].header
        model_image     = img[2].data
        model_header    = img[2].header
        residual_image  = img[3].data
        residual_header = img[3].header

        sb_original_image = create_sb_image(original_image,zp,pix_scale)
        sb_model_image    = create_sb_image(model_image,   zp,pix_scale)
        sb_residual_image = create_sb_image(residual_image,zp,pix_scale)

        fits.writeto("original_"+type_of_analysis+"_sb.fits", sb_original_image,original_header,clobber=True)
        fits.writeto("model_"+type_of_analysis+"_sb.fits",       sb_model_image,original_header,clobber=True) #because original header has all the interesting keywords for plotting
        fits.writeto("residuals_"+type_of_analysis+"_sb.fits",sb_residual_image,original_header,clobber=True) #because original header has all the interesting keywords for plotting
        fits.writeto("residuals_"+type_of_analysis+"_no_sb.fits",residual_image,original_header,clobber=True) #because original header has all the interesting keywords for plotting

        img_original_sb = fits.open("original_"+type_of_analysis+"_sb.fits")
        img_model_sb    = fits.open("model_"+type_of_analysis+"_sb.fits")
        img_residual_sb = fits.open("residuals_"+type_of_analysis+"_sb.fits")
        img_residual    = fits.open("residuals_"+type_of_analysis+"_no_sb.fits")

        if type_of_analysis == "single":
            re_pix,re_kpc,nn,ar,pa = get_struct_param_best_fit(gal,band,type_of_analysis)
        else:
            re1_pix,re1_kpc,nn1,ar1,pa1,re2_pix,re2_kpc,nn2,ar2,pa2 = get_struct_param_best_fit(gal,band,type_of_analysis)

        subfig = aplpy.FITSFigure(img_original_sb, figure = fig, subplot=plotting_pos[cont], origin='lower')
        y0,x0 = finding_coo_brightest_SB(img_original_sb)
        max_SB = img_original_sb[0].data[x0,y0]
        min_SB = noise_level
        plotting_options()
        subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
        if flag_final_sample == True:
            if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
        plotting_data_original(gal,band,filename,zz[cont_gal],logmass[cont_gal])
        plotting_scalebar()
        cont = cont+1

        subfig = aplpy.FITSFigure(img_model_sb, figure = fig, subplot=plotting_pos[cont], origin='lower')
        if type_of_analysis == "single":
            printing_effective_radius(img_model_sb,pix_scale,re_pix*pix_scale,ar,pa) #only in the original image, as it is the only one with "good" header
            plotting_data_single_model(re_kpc,re_pix,nn,ar)
        else:
            plotting_data_double_model(re1_kpc,re1_pix,nn1,ar1,re2_kpc,re2_pix,nn2,ar2)
        plotting_options()
        subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
        if flag_final_sample == True:
            if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
        cont = cont+1

        subfig = aplpy.FITSFigure(img_residual_sb, figure = fig, subplot=plotting_pos[cont], origin='lower')
        plotting_options()
        subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
        if flag_final_sample == True:
            if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
        #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        box_position = [0.275,0.8-(cont_band*0.25),0.2,0.025] #[xmin,ymin,dx,dy]
        plotting_colorbar(box_position,max_SB,min_SB) #only for the last image, create the colorbar
        #=============================
        #the colorbar changes the subimage layout, I need to print it again
        subfig = aplpy.FITSFigure(img_residual_sb, figure = fig, subplot=plotting_pos[cont], origin='lower')
        plotting_options()
        subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
        if flag_final_sample == True:
            if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
        #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        #=============================
        cont = cont+1

        subfig = aplpy.FITSFigure(img_residual, figure = fig, subplot=plotting_pos[cont], origin='lower')
        plotting_options()
        subfig.show_colorscale(cmap='cubehelix')
        if flag_final_sample == True:
            if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
        #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        box_position = [0.775,0.30,0.2,0.025] #[xmin,ymin,dx,dy]
        if band == bands[-1]: plotting_colorbar(box_position) #only for the last image, create the colorbar
        #=============================
        #the colorbar changes the subimage layout, I need to print it again
        if band == bands[-1]:
            subfig = aplpy.FITSFigure(img_residual, figure = fig, subplot=plotting_pos[cont], origin='lower')
            plotting_options()
            subfig.show_colorscale(cmap='cubehelix')
            if flag_final_sample == True:
                if vector_good_fit[cont_band] == False: subfig.show_rectangles(ra[filter_big_cat], dec[filter_big_cat],1,1,facecolor="red",alpha=0.5) #if the analysis was wrong, it is highlighted in red color (as the stop in a semaphore)
            #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        #=============================
        cont = cont+1
    
    fig.canvas.draw()
    plt.savefig("./individual_galaxy_SB_profiles/"+str(gal)+"_galfit_fits_"+type_of_analysis+".pdf", format='pdf', dpi=100)
    plt.clf()
    plt.close(fig)

    os.remove("original_"+type_of_analysis+"_sb.fits")
    os.remove("model_"+type_of_analysis+"_sb.fits")
    os.remove("residuals_"+type_of_analysis+"_sb.fits")
    os.remove("residuals_"+type_of_analysis+"_no_sb.fits")
