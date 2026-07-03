# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 10:24:28 2023

@author: DELINTE Nicolas
"""

import os
import sys
import json
import warnings
import numpy as np
import nibabel as nib
from scipy.linalg import polar
from unravel.utils import tensor_to_DTI, get_streamline_density, tract_to_ROI
from unravel.core import (get_fixel_weight, get_microstructure_map,
                          get_weighted_mean, tensor_to_peak)
from unravel.analysis import get_metric_along_trajectory
from unravel.stream import extract_nodes, get_roi_sections_from_nodes

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from dipy.io.streamline import load_tractogram


def to_float64(val):
    """
    Used if *val* is an instance of numpy.float32.
    """

    return np.float64(val)


def create_tensor_metrics(path: str):
    '''


    Parameters
    ----------
    path : str
        Ex: '/.../diamond/subjectName'

    Returns
    -------
    None.

    '''

    metric = {}
    for tensor in ['t0', 't1']:

        img = nib.load(path + '_diamond_' + tensor + '.nii.gz')
        t = img.get_fdata()

        FA, AD, RD, MD = tensor_to_DTI(t)

        metric['FA_'+tensor] = FA
        metric['MD_'+tensor] = MD
        metric['AD_'+tensor] = AD
        metric['RD_'+tensor] = RD

        for m in metric:
            out = nib.Nifti1Image(metric[m].real, img.affine)
            out.header.get_xyzt_units()
            out.to_filename(path + '_diamond_' + m + '.nii.gz')

    np.seterr(divide='ignore', invalid='ignore')

    metric_list = ['FA', 'MD', 'AD', 'RD']
    metrics = {}

    fracs = nib.load(path+'_diamond_fractions.nii.gz').get_fdata()

    for comp in range(3):

        metrics['fractions_t'+str(comp)] = fracs[:, :, :, 0, comp]

    for i in metric_list:
        metric['w'+i] = (metric[i+'_t0']*metrics['fractions_t0'] +
                         metric[i+'_t1']*metrics['fractions_t1'])/(
            metrics['fractions_t0']+metrics['fractions_t1'])
        metric['w'+i][np.isnan(metric['w'+i])] = 0

    np.seterr(divide='warn', invalid='warn')

    for m in list(metric.keys()):
        print(m)
        if 'w' in m:
            out = nib.Nifti1Image(metric[m], img.affine, img.header)
            out.header.get_xyzt_units()
            out.to_filename(path+'_diamond_'+m+'.nii.gz')


def get_mean_tracts(trk_file: str, micro_path: str, trajectory: bool = False,
                    method: str = 'ang', weighting: str = 'tsl'):
    '''
    Return means for all metrics for a single patient using UNRAVEL

    Parameters
    ----------
    trk_file : str
        DESCRIPTION.
    micro_path : str
        Patient specific path to microstructure folder

    Returns
    -------
    mean : TYPE
        DESCRIPTION.
    dev : TYPE
        DESCRIPTION.

    '''

    trk = load_tractogram(trk_file, 'same')
    trk.to_vox()
    trk.to_corner()

    subject = micro_path.split('/')[-4]

    mean_dic = {}
    dev_dic = {}

    # Streamlines properties ---------------------

    if trajectory:
        print(trk_file)
        point_array = extract_nodes(trk_file)
        roi_sections = get_roi_sections_from_nodes(trk_file, point_array)
        mean_dic['voxel_count'] = np.unique(roi_sections,
                                            return_counts=True)[1].tolist()
        dev_dic['voxel_count'] = 0
    else:
        mean_dic['stream_count'] = len(trk.streamlines._offsets)
        dev_dic['stream_count'] = 0
        mean_dic['voxel_count'] = np.sum(tract_to_ROI(trk_file))
        dev_dic['voxel_count'] = 0

    # DTI ----------------------------------

    if os.path.isdir(micro_path + 'dti/'):

        metric_list = ['FA', 'AD', 'MD', 'RD']

        for m in metric_list:

            map_file = micro_path + 'dti/' + subject + '_' + m + '.nii.gz'

            metric_maps = nib.load(map_file).get_fdata()[..., np.newaxis]

            fixel_weights = get_streamline_density(trk)[..., np.newaxis]

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = metric_maps[..., -1]
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m+'_DTI'] = mean
            dev_dic[m+'_DTI'] = dev

    # NODDI -------------------------------------

    if os.path.isdir(micro_path + 'noddi/'):

        metric_list = ['fiso', 'fintra', 'fextra', 'odi']

        for m in metric_list:

            map_file = micro_path + 'noddi/' + subject + '_noddi_' + m + '.nii.gz'

            metric_maps = nib.load(map_file).get_fdata()[..., np.newaxis]

            fixel_weights = get_streamline_density(trk)[..., np.newaxis]

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = metric_maps[..., -1]
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    # Diamond (wM) --------------------------------

    if os.path.isdir(micro_path + 'diamond/') and os.path.isfile(micro_path + 'diamond/' + subject + '_diamond_wFA.nii.gz'):

        if not os.path.isfile(micro_path + 'diamond/' + subject
                              + '_diamond_FA_t0.nii.gz'):

            create_tensor_metrics(micro_path + 'diamond/' + subject)

        metric_list = ['wFA', 'wAD', 'wMD', 'wRD']

        for m in metric_list:

            map_file = micro_path + 'diamond/' + subject + '_diamond_' + m + '.nii.gz'

            metric_maps = nib.load(map_file).get_fdata()[..., np.newaxis]

            fixel_weights = get_streamline_density(trk)[..., np.newaxis]

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = metric_maps[..., -1]
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    # Diamond ------------------------------

    if os.path.isdir(micro_path + 'diamond/'):

        tensor_files = [micro_path + 'diamond/' + subject + '_diamond_t0.nii.gz',
                        micro_path + 'diamond/' + subject + '_diamond_t1.nii.gz']

        peaks = np.stack((tensor_to_peak(nib.load(tensor_files[0]).get_fdata()),
                          tensor_to_peak(nib.load(tensor_files[1]).get_fdata())),
                         axis=4)

        fixel_weights = get_fixel_weight(trk, peaks, method=method)

        fracs = nib.load(micro_path + 'diamond/' + subject
                         + '_diamond_fractions.nii.gz').get_fdata()
        metric_maps = np.stack(
            (fracs[:, :, :, 0, 0], fracs[:, :, :, 0, 1]), axis=3)

        if trajectory:
            mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                    roi_sections,
                                                    weighting=weighting)
            mean = mean.tolist()
            dev = dev.tolist()
        else:
            microstructure_map = get_microstructure_map(
                fixel_weights, metric_maps)
            mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                          weighting=weighting)

        mean_dic['frac_dmd'] = mean
        dev_dic['frac_dmd'] = dev

        metric_list = ['FA', 'MD', 'RD', 'AD']

        if not os.path.isfile(micro_path + 'diamond/' + subject
                              + '_diamond_FA_t0.nii.gz'):

            create_tensor_metrics(micro_path + 'diamond/' + subject)

        for m in metric_list:

            map_files = [micro_path + 'diamond/' + subject + '_diamond_' + m
                         + '_t0.nii.gz', micro_path + 'diamond/' + subject
                         + '_diamond_' + m + '_t1.nii.gz']

            metric_maps = np.stack((nib.load(map_files[0]).get_fdata(),
                                    nib.load(map_files[1]).get_fdata()),
                                   axis=3)

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = get_microstructure_map(fixel_weights,
                                                            metric_maps)
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    # MF (wM-ish) ------------------------------

    if os.path.isdir(micro_path + 'mf/'):

        metric_list = ['frac_csf', 'fvf_tot']

        for m in metric_list:

            map_file = micro_path + 'mf/' + subject + '_mf_' + m + '.nii.gz'

            metric_maps = nib.load(map_file).get_fdata()[..., np.newaxis]

            fixel_weights = get_streamline_density(trk)[..., np.newaxis]

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = metric_maps[..., -1]
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    # Microstructure fingerprinting --------

    if os.path.isdir(micro_path + 'mf/'):

        tensor_files = [micro_path + 'mf/' + subject + '_mf_peak_f0.nii.gz',
                        micro_path + 'mf/' + subject + '_mf_peak_f1.nii.gz']

        peaks = np.stack((nib.load(tensor_files[0]).get_fdata(),
                          nib.load(tensor_files[1]).get_fdata()),
                         axis=4)

        fixel_weights = get_fixel_weight(trk, peaks, method=method)

        metric_list = ['fvf', 'frac']

        for m in metric_list:

            map_files = [micro_path + 'mf/' + subject + '_mf_' + m + '_f0.nii.gz',
                         micro_path + 'mf/' + subject + '_mf_' + m + '_f1.nii.gz']

            metric_maps = np.stack((nib.load(map_files[0]).get_fdata(),
                                    nib.load(map_files[1]).get_fdata()),
                                   axis=3)

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights, metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = get_microstructure_map(fixel_weights,
                                                            metric_maps)
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    # AFD ----------------------------------------------------------------

    if os.path.isdir(micro_path + 'afd/'):

        peaks_img = nib.load(micro_path + 'afd/'+subject+'_peaks.nii.gz')
        peaks = peaks_img.get_fdata()

        u, _ = polar(peaks_img.affine[0:3, 0:3])
        K = int(peaks.shape[-1]/3)

        peaks_unravel = np.zeros((peaks.shape[:-1])+(3, K,))
        for k in range(K):
            peaks_unravel[..., :, k] = peaks[..., 3*k:3*k+3] @ u

        fixel_weights = get_fixel_weight(trk, peaks_unravel, method=method)

        metric_list = ['afd', 'disp']

        for m in metric_list:

            metric_maps = nib.load(micro_path + 'afd/' +
                                   subject+'_'+m+'_voxel.nii.gz').get_fdata()

            if trajectory:
                mean, dev = get_metric_along_trajectory(fixel_weights,
                                                        metric_maps,
                                                        roi_sections,
                                                        weighting=weighting)
                mean = mean.tolist()
                dev = dev.tolist()
            else:
                microstructure_map = get_microstructure_map(fixel_weights,
                                                            metric_maps)
                mean, dev = get_weighted_mean(microstructure_map, fixel_weights,
                                              weighting=weighting)

            mean_dic[m] = mean
            dev_dic[m] = dev

    return mean_dic, dev_dic


def get_mean_tracts_study(root: str, region_list: list,
                          output_path: str, subj_list: list = None,
                          trajectory: bool = False, method: str = 'ang',
                          weighting: str = 'tsl'):
    '''


    Parameters
    ----------
    root : str
        DESCRIPTION.
    selected_edges_path : str
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    subjects_list = root + 'subjects/subj_list.json'

    if subj_list is None:
        with open(subjects_list, 'r') as read_file:
            subj_list = json.load(read_file)

    dic_tot = {}
    dic_tot['Mean'] = {}
    dic_tot['Dev'] = {}

    for sub in subj_list:

        micro_path = root + 'subjects/' + sub + '/dMRI/microstructure/'
        tract_path = root + 'subjects/' + sub + '/dMRI/tractography/tois/'

        dic_tot['Mean'][sub] = {}
        dic_tot['Dev'][sub] = {}

        eddy_qc_file = (root+'subjects/'+sub+'/dMRI/preproc/eddy/'+sub
                        + '_eddy_corr.qc/qc.json')

        if os.path.isfile(eddy_qc_file):
            print('File found for movement metrics')

            qc = json.load(open(eddy_qc_file))

            move = qc['qc_mot_rel']
            snr = qc['qc_cnr_avg'][0]

        for roi in region_list:

            try:
                trk_file = (tract_path + sub + '_' + roi + '.trk')

                mean_dic, dev_dic = get_mean_tracts(trk_file, micro_path,
                                                    trajectory=trajectory,
                                                    method=method,
                                                    weighting=weighting)

            except FileNotFoundError:
                print('.trk file or metrics not found for region ' + str(roi)
                      + ' in patient ' + sub + ' at '+trk_file)
                continue
            except IndexError:
                print('IndexError with subject ' + sub)
                continue
            except ValueError:
                print('Trajectory: Insufficient streamlines for region '
                      + str(roi) + ' in patient ' + sub)

            if os.path.isfile(eddy_qc_file):
                mean_dic['movement'] = move
                mean_dic['snr'] = snr

            dic_tot['Mean'][sub][roi] = mean_dic
            dic_tot['Dev'][sub][roi] = dev_dic

    if trajectory:
        dic_file = (output_path + 'unravel_trajectory_'+sub
                    + '_'+method+'_'+weighting+'.json')
    else:
        dic_file = output_path + 'unravel_'+sub+'_'+method+'_'+weighting+'.json'
    json.dump(dic_tot.copy(), open(dic_file, 'w'),
              default=to_float64, indent=4)


def merge_json_dics(folder_path, trajectory: bool = False,
                    method: str = 'ang', weighting: str = 'tsl'):

    if trajectory:
        dic_file = output_path + 'unravel_trajectory'+'_'+method+'_'+weighting+'.json'
    else:
        dic_file = output_path + 'unravel_mean'+'_'+method+'_'+weighting+'.json'
    if os.path.isfile(dic_file):
        os.system('rm '+dic_file)

    merged_dict = {}
    merged_dict['Mean'] = {}
    merged_dict['Dev'] = {}

    # Iterate through files in the directory
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        # Check if the file is a JSON file
        if filename.endswith('.json'):
            if weighting in filename and method in filename:
                if trajectory:
                    if 'trajectory' not in filename:
                        continue
                else:
                    if 'trajectory' in filename:
                        continue
                # try:
                with open(file_path, 'r') as file:
                    data = json.load(file)
                    merged_dict['Mean'].update(data['Mean'])
                    merged_dict['Dev'].update(data['Dev'])
                # except Exception as e:
                #     print(f"Error reading file '{file_path}': {e}")

    json.dump(merged_dict, open(dic_file, 'w'), default=to_float64, indent=4)

    return merged_dict


if __name__ == '__main__':

    root = sys.argv[1]

    method = 'ang'
    weighting = 'tsl'
    merge = (sys.argv[2] == "merge")

    region_list = ['cc_ant_midbody', 'cc_post_midbody', 'cc_genu', 'cc_isthmus',
                   'cc_splenium', 'uf_left', 'atr_left', 'fornix_left', 'uf_right',
                   'atr_right', 'fornix_right', 'cingulum_left', 'cingulum_right',
                   'ci_left', 'ci_right']

    output_path = sys.argv[1]+'out/'

    os.makedirs(output_path, exist_ok=True)

    if not merge:
        get_mean_tracts_study(root, region_list, output_path,
                              subj_list=[sys.argv[2]],
                              trajectory=False, method=method,
                              weighting=weighting)
        get_mean_tracts_study(root, region_list, output_path,
                              subj_list=[sys.argv[2]],
                              trajectory=True, method=method,
                              weighting=weighting)

    else:
        # Merge dictionaries from JSON files in the folder
        merged_dictionary = merge_json_dics(output_path, trajectory=False,
                                            method=method, weighting=weighting)
        merged_dictionary = merge_json_dics(output_path, trajectory=True,
                                            method=method, weighting=weighting)
