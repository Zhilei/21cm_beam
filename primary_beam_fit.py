import numpy as np
import matplotlib.pyplot as plt
from pyuvdata import UVData
import healpy as hp
from glob import glob
import hera_cal
import copy

from direct_optimal_mapping import data_conditioning, pixel_selection
from direct_optimal_mapping import optimal_mapping_radec_grid

import beam_functions as func
from scipy import optimize

import hera_cal
import pandas as pd
from astropy import constants

ra_ctr_deg = 30
ra_rng_deg = 20

bl_max = 300 # m

src_ra_deg = 30.3
src_dec_deg = -30.8

ipol = -5

data_files = np.array(sorted(glob('/nfs/esc/hera/H6C/IDR2/2459873/zen.2459873.*.sum.omni_vis.uvh5')))

jds = np.array([float(data_file.split('/')[-1][4:17]) for data_file in data_files])
# print(data_files)

lsts = hera_cal.utils.JD2LST(jds)

idx = np.where((lsts > np.radians(ra_ctr_deg - ra_rng_deg/2.)) & (lsts < np.radians(ra_ctr_deg + ra_rng_deg/2.)))[0]
print('%d files selected.'%len(idx))



for ifreq in range(1400, 100, -100):
    z_pointing = []
    fit_result = []
    for data_file in data_files[idx]:
        print('File:', data_file)
        uv = UVData()
        uv.read(data_file)
        time_arr = np.unique(uv.time_array)
        uv.unphase_to_drift()
        uv.phase_type = 'drift'

        for time_idx in range(len(time_arr)):

            uv_sel = uv.select(times=time_arr[time_idx], inplace=False, keep_all_metadata=False)
            print('################\nSelected data shape:', uv_sel.data_array.shape)
            lst_ze = np.unique(uv_sel.lst_array)[0]
            print('Zenith pointing:', np.degrees(lst_ze))

            dc = data_conditioning.DataConditioning(uv_sel, ifreq, ipol)
            print('Frequency: %.2fMHz'%(dc.uv_1d.freq_array[0]/1e6))
            dc.noise_calc()
            dc.rm_flag()
            dc.redundant_avg()
            print('Data array shape after data conditioning:', dc.uv_1d.data_array.shape)
            print('Noise array shape after data conditioning:', dc.uvn.data_array.shape)
            syn_beam = np.degrees(constants.c.value/dc.uv_1d.freq_array[0][0]/bl_max)

            # Pixels Setup
            ra_ctr_deg = src_ra_deg 
            ra_rng_deg = syn_beam * 5
            n_ra = 40
            dec_ctr_deg = src_dec_deg
            dec_rng_deg = syn_beam * 5 
            n_dec = 40

            sky_px = optimal_mapping_radec_grid.SkyPx()
            px_dic = sky_px.calc_radec_pix(ra_ctr_deg, ra_rng_deg, n_ra, dec_ctr_deg, dec_rng_deg, n_dec)
            shape = px_dic['sa_sr'].shape

            # Optimal Mapping & A matrix
            opt_map = optimal_mapping_radec_grid.OptMapping(dc.uv_1d, px_dic, epoch='Current')
            opt_map.set_a_mat()
            inv_noise_mat = opt_map.set_inv_noise_mat(dc.uvn, norm=True)

            print('a_mat shape is:', opt_map.a_mat.shape)

            map_vis = np.matmul(np.conjugate(opt_map.a_mat.T), 
                                  np.matmul(opt_map.inv_noise_mat, 
                                            np.matrix(opt_map.data)))

            # map_vis = np.real(map_vis)
            map_vis = np.array(map_vis.real.reshape(px_dic['px_id'].shape))#/px_dic['sa_sr']

            weight = np.matmul(np.conjugate(opt_map.beam_mat.T), np.diag(opt_map.inv_noise_mat)).reshape(px_dic['px_id'].shape)
            weight_sq = np.matmul(np.conjugate((opt_map.beam_mat**2).T), np.diag(opt_map.inv_noise_mat)).reshape(px_dic['px_id'].shape)

            map_norm = map_vis/weight

            xieta = np.radians(list(zip(px_dic['ra_deg'].flatten(), px_dic['dec_deg'].flatten()))).T

            # initial guesses then fitting
            p0 = [15, 
                np.radians(src_ra_deg), 
                np.radians(src_dec_deg), 
                np.radians(0.5), 
                np.radians(0.5), 
                np.radians(0)]

            bounds = ([0, np.radians(src_ra_deg-1), np.radians(src_dec_deg-1), np.radians(0.5 * syn_beam), np.radians(0.5 * syn_beam), -np.pi],
                      [1000, np.radians(src_ra_deg+1), np.radians(src_dec_deg+1), np.radians(2.0 * syn_beam), np.radians(2.0 * syn_beam), np.pi],)

            try:
                popt, pcov = optimize.curve_fit(func.gaussian2d, xieta, map_norm.flatten(), p0=p0, bounds=bounds)
            except:
                print('Fitting Failed.')
                continue

            print(popt)
            z_pointing.append(lst_ze)
            fit_result.append(popt)

    z_pointing = np.array(z_pointing)
    fit_result = np.array(fit_result)
    result = np.concatenate((z_pointing[:, np.newaxis], fit_result), axis=1)
    df = pd.DataFrame(data=result, columns=['Pointing', 'Amp', 'x0', 'y0', 'fwhm_major', 'fwhm_minor', 'fwhm_theta'])
    df.to_csv('primary_beam_%6.2fMHz.csv'%(uv.freq_array[0, ifreq]/1e6), index=False)