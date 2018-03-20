import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import os
import pdb


def bool2int(ele):
    if ele == True:
        return 1
    else:
        return 0


bands=["g","r","i","z"]
type_of_analysis = "single"
min_number_good_and_small = 2
galaxies_to_be_excluded = ["373300","537226"] #large galaxies with compact nucleus/bulge
#selection criteria
value_sel_flag_obj_detected = 0
value_max_size = 2.
value_min_size = 0.01
value_sel_re_exist = 1
value_sel_max_nn = 10
value_sel_min_nn = 0.1
value_sel_nn_exist = 1
value_sel_ar_exist = 1


tt_mass = Table.read("./summary_masses_for_sample262objs_Ferreras_update.cat",
                     names = ["cataid","z_tonry","nbands","s2n","logmstar","err_logmstar","logage","err_logage","metal","err_metal","url","logmstar_Ferreras"],
                     format="ascii.commented_header")

#reading the catalogs (individually to make things faster, in comparison to use a function)
ttg_single = Table.read("./analysis_filter_g/gama_kids_g_single_sersic_fits.cat",\
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],\
             format="ascii.commented_header")

ttr_single = Table.read("./analysis_filter_r/gama_kids_r_single_sersic_fits.cat",
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],\
             format="ascii.commented_header")

tti_single = Table.read("./analysis_filter_i/gama_kids_i_single_sersic_fits.cat",
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],\
             format="ascii.commented_header")
     
ttz_single = Table.read("./analysis_filter_z/viking_z_single_sersic_fits.cat",
             names = ["gal_id","zz","chi2","chi2nu","file_to_restart","x_center","y_center",\
                          "mag","mag_error","mag_exist",\
                          "re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err",\
                          "n","n_error","n_exist",\
                          "ar","ar_error","ar_exist",\
                          "pa","pa_error","pa_exist",\
                          "index_of_best_fit","flag_obj_detected"],\
             format="ascii.commented_header")
     
filt_obj_det_g_single  = ttg_single["flag_obj_detected"] == value_sel_flag_obj_detected
filt_max_size_g_single = ttg_single["re_kpc"] > value_max_size
filt_min_size_g_single = ttg_single["re_kpc"] < value_min_size
filt_re_exist_g_single = ttg_single["re_exist"] == value_sel_re_exist
filt_max_nn_g_single   = ttg_single["n"] == value_sel_max_nn
filt_min_nn_g_single   = ttg_single["n"] == value_sel_min_nn
filt_nn_exist_g_single = ttg_single["n_exist"] == value_sel_nn_exist
filt_ar_exist_g_single = ttg_single["ar_exist"] == value_sel_ar_exist

filt_obj_det_r_single  = ttr_single["flag_obj_detected"] == value_sel_flag_obj_detected
filt_max_size_r_single = ttr_single["re_kpc"] > value_max_size
filt_min_size_r_single = ttr_single["re_kpc"] < value_min_size
filt_re_exist_r_single = ttr_single["re_exist"] == value_sel_re_exist
filt_max_nn_r_single   = ttr_single["n"] == value_sel_max_nn
filt_min_nn_r_single   = ttr_single["n"] == value_sel_min_nn
filt_nn_exist_r_single = ttr_single["n_exist"] == value_sel_nn_exist
filt_ar_exist_r_single = ttr_single["ar_exist"] == value_sel_ar_exist

filt_obj_det_i_single  = tti_single["flag_obj_detected"] == value_sel_flag_obj_detected
filt_max_size_i_single = tti_single["re_kpc"] > value_max_size
filt_min_size_i_single = tti_single["re_kpc"] < value_min_size
filt_re_exist_i_single = tti_single["re_exist"] == value_sel_re_exist
filt_max_nn_i_single   = tti_single["n"] == value_sel_max_nn
filt_min_nn_i_single   = tti_single["n"] == value_sel_min_nn
filt_nn_exist_i_single = tti_single["n_exist"] == value_sel_nn_exist
filt_ar_exist_i_single = tti_single["ar_exist"] == value_sel_ar_exist

filt_obj_det_z_single  = ttz_single["flag_obj_detected"] == value_sel_flag_obj_detected
filt_max_size_z_single = ttz_single["re_kpc"] > value_max_size
filt_min_size_z_single = ttz_single["re_kpc"] < value_min_size
filt_re_exist_z_single = ttz_single["re_exist"] == value_sel_re_exist
filt_max_nn_z_single   = ttz_single["n"] == value_sel_max_nn
filt_min_nn_z_single   = ttz_single["n"] == value_sel_min_nn
filt_nn_exist_z_single = ttz_single["n_exist"] == value_sel_nn_exist
filt_ar_exist_z_single = ttz_single["ar_exist"] == value_sel_ar_exist

""" before 15-Jan-2018
#the next np.logical_or.reduce is the way to generalize the logical_or for several arrays
total_filt_bad_g  = np.logical_or.reduce((filt_obj_det_g_single,filt_max_size_g_single,filt_min_size_g_single,filt_re_exist_g_single,filt_max_nn_g_single,filt_min_nn_g_single,filt_nn_exist_g_single,filt_ar_exist_g_single),dtype=bool)
total_filt_bad_g_ids  = np.array(ttg_single["gal_id"])[total_filt_bad_g]
total_filt_good_g = np.logical_not(total_filt_bad_g)
total_filt_good_g_ids = np.array(ttg_single["gal_id"])[total_filt_good_g]

total_filt_bad_r  = np.logical_or.reduce((filt_obj_det_r_single,filt_max_size_r_single,filt_min_size_r_single,filt_re_exist_r_single,filt_max_nn_r_single,filt_min_nn_r_single,filt_nn_exist_r_single,filt_ar_exist_r_single),dtype=bool)
total_filt_bad_r_ids  = np.array(ttr_single["gal_id"])[total_filt_bad_r]
total_filt_good_r = np.logical_not(total_filt_bad_r)
total_filt_good_r_ids = np.array(ttr_single["gal_id"])[total_filt_good_r]

total_filt_bad_i  = np.logical_or.reduce((filt_obj_det_i_single,filt_max_size_i_single,filt_min_size_i_single,filt_re_exist_i_single,filt_max_nn_i_single,filt_min_nn_i_single,filt_nn_exist_i_single,filt_ar_exist_i_single),dtype=bool)
total_filt_bad_i_ids  = np.array(tti_single["gal_id"])[total_filt_bad_i]
total_filt_good_i = np.logical_not(total_filt_bad_i)
total_filt_good_i_ids = np.array(tti_single["gal_id"])[total_filt_good_i]

total_filt_bad_z  = np.logical_or.reduce((filt_obj_det_z_single,filt_max_size_z_single,filt_min_size_z_single,filt_re_exist_z_single,filt_max_nn_z_single,filt_min_nn_z_single,filt_nn_exist_z_single,filt_ar_exist_z_single),dtype=bool)
total_filt_bad_z_ids  = np.array(ttz_single["gal_id"])[total_filt_bad_z]
total_filt_good_z = np.logical_not(total_filt_bad_z)
total_filt_good_z_ids = np.array(ttz_single["gal_id"])[total_filt_good_z]

print(total_filt_good_g_ids.size,total_filt_good_r_ids.size,total_filt_good_i_ids.size,total_filt_good_z_ids.size) #how many valid galaxies per filter

total_filt = np.logical_or.reduce((total_filt_good_g,total_filt_good_r,total_filt_good_i,total_filt_good_z),dtype=bool)
print(np.count_nonzero(total_filt==True)) #to see how many galaxies in total (redundant)
"""

#good and bad fits
total_filt_bad_g  = np.logical_or.reduce((filt_obj_det_g_single,filt_min_size_g_single,filt_re_exist_g_single,filt_max_nn_g_single,filt_min_nn_g_single,filt_nn_exist_g_single,filt_ar_exist_g_single),dtype=bool)
total_filt_good_g = np.logical_not(total_filt_bad_g)
total_filt_bad_r  = np.logical_or.reduce((filt_obj_det_r_single,filt_min_size_r_single,filt_re_exist_r_single,filt_max_nn_r_single,filt_min_nn_r_single,filt_nn_exist_r_single,filt_ar_exist_r_single),dtype=bool)
total_filt_good_r = np.logical_not(total_filt_bad_r)
total_filt_bad_i  = np.logical_or.reduce((filt_obj_det_i_single,filt_min_size_i_single,filt_re_exist_i_single,filt_max_nn_i_single,filt_min_nn_i_single,filt_nn_exist_i_single,filt_ar_exist_i_single),dtype=bool)
total_filt_good_i = np.logical_not(total_filt_bad_i)
total_filt_bad_z  = np.logical_or.reduce((filt_obj_det_z_single,filt_min_size_z_single,filt_re_exist_z_single,filt_max_nn_z_single,filt_min_nn_z_single,filt_nn_exist_z_single,filt_ar_exist_z_single),dtype=bool)
total_filt_good_z = np.logical_not(total_filt_bad_z)

#good and small fits
total_filt_bad_small_g  = np.logical_or.reduce((filt_obj_det_g_single,filt_max_size_g_single,filt_min_size_g_single,filt_re_exist_g_single,filt_max_nn_g_single,filt_min_nn_g_single,filt_nn_exist_g_single,filt_ar_exist_g_single),dtype=bool)
total_filt_good_small_g = np.logical_not(total_filt_bad_small_g)
total_filt_bad_small_r  = np.logical_or.reduce((filt_obj_det_r_single,filt_max_size_r_single,filt_min_size_r_single,filt_re_exist_r_single,filt_max_nn_r_single,filt_min_nn_r_single,filt_nn_exist_r_single,filt_ar_exist_r_single),dtype=bool)
total_filt_good_small_r = np.logical_not(total_filt_bad_small_r)
total_filt_bad_small_i  = np.logical_or.reduce((filt_obj_det_i_single,filt_max_size_i_single,filt_min_size_i_single,filt_re_exist_i_single,filt_max_nn_i_single,filt_min_nn_i_single,filt_nn_exist_i_single,filt_ar_exist_i_single),dtype=bool)
total_filt_good_small_i = np.logical_not(total_filt_bad_small_i)
total_filt_bad_small_z  = np.logical_or.reduce((filt_obj_det_z_single,filt_max_size_z_single,filt_min_size_z_single,filt_re_exist_z_single,filt_max_nn_z_single,filt_min_nn_z_single,filt_nn_exist_z_single,filt_ar_exist_z_single),dtype=bool)
total_filt_good_small_z = np.logical_not(total_filt_bad_small_z)

#the condition to pertain to the final sample is that the galaxy has good and small analysis in at least one filter
#total_filt = np.logical_or.reduce((total_filt_good_small_g,total_filt_good_small_r,total_filt_good_small_i,total_filt_good_small_z),dtype=bool)
#print(np.count_nonzero(total_filt==True)) #to see how many galaxies in total (redundant)

#the condition to pertain to the final sample is that the galaxy has good and small analysis in at least min_number_good_and_small filters
how_many_good = np.zeros(len(ttg_single["gal_id"]))
for ii in range(len(ttg_single["gal_id"])):
    how_many_good[ii] = bool2int(total_filt_good_small_g[ii]) + bool2int(total_filt_good_small_r[ii]) + bool2int(total_filt_good_small_i[ii]) + bool2int(total_filt_good_small_z[ii])
    if str(ttg_single["gal_id"][ii]) in galaxies_to_be_excluded: how_many_good[ii] = 0
    print(str(ttg_single["gal_id"][ii]),how_many_good[ii])
total_filt = how_many_good >= min_number_good_and_small
print(np.count_nonzero(total_filt==True)) #to see how many galaxies in total

#rejecting objects that are not massive enough
mass = 10.**(tt_mass["logmstar_Ferreras"])
total_filt_mass = mass > 8e10
total_filt = np.logical_and(total_filt,total_filt_mass)

#borderline objects - uncomment if necessary
"""
filt1 = mass > 6e10
filt2 = mass < 8e10
total_filt = np.logical_and.reduce((total_filt,filt1,filt2),dtype=bool)
"""

galaxy_name = np.array(ttg_single["gal_id"])[total_filt]
galaxy_name = galaxy_name.astype(str)
galaxy_zz   = np.array(ttg_single["zz"])[total_filt]
#I order the galaxies according to their redshift
zz_ordered_indices = np.argsort(galaxy_zz)
galaxy_name = galaxy_name[zz_ordered_indices]

command  = "pdfunite "
command2 = "pdfunite "

for gal in galaxy_name:
    command  = command + "./individual_galaxy_SB_profiles/"+ gal + "_galfit_fits_"+type_of_analysis+".pdf "
    command2 = command2 + "./color_images/"+gal+"_rgb.pdf "
    
print(command  + "all_gama_"+type_of_analysis+"_good.pdf")
print(" ")
print(command2 + "all_gama_"+type_of_analysis+"_good_rgb.pdf")
print(" ")
os.system(command  + "all_gama_"+type_of_analysis+"_good.pdf")
os.system(command2 + "all_gama_"+type_of_analysis+"_good_rgb.pdf")

#file_to_write = open("list_good_"+type_of_analysis+"_gals_borderline.txt","w")
file_to_write = open("list_good_"+type_of_analysis+"_gals.txt","w")
for gal in galaxy_name:
    file_to_write.write(gal+"\n")
file_to_write.close()
