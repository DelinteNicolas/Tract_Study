# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:48:12 2023

@author: N. Delinte

Code used to track specific white matter bundles.
"""

import os
import sys
import json
import numpy as np
import nibabel as nib
from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from regis.core import find_transform, apply_transform
from unravel.stream import extract_nodes, remove_outlier_streamlines


def register_rois_to_subj(fa_path: str, rois_path: str, mni_fa_path: str,
                          mni_T1_path: str, T1_path: str, output_path: str,
                          static_mask_path: str, T1_based_regis: bool = False,
                          binary_thresh: float = 0.25):
    '''
    Two-step registration to obtain label in the diffusion space.

    Parameters
    ----------
    fa_path : str
        Path to FA file.
    roi_path : str
        DESCRIPTION.
    mni_fa_path : str
        DESCRIPTION.
    output_path : str
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    if T1_based_regis:

        map_mni_to_subj_T1 = find_transform(mni_T1_path, T1_path)
        map_T1_to_subj = find_transform(T1_path, fa_path, only_affine=True)
        static = nib.load(fa_path)

        for roi_path in os.listdir(rois_path):

            if '.nii.gz' not in roi_path:
                continue

            output_filename = output_path+roi_path.split('/')[-1]

            T1_step = apply_transform(rois_path+roi_path, map_mni_to_subj_T1,
                                      binary=True, binary_thresh=binary_thresh)
            transformed = map_T1_to_subj.transform(T1_step)
            transformed = np.where(transformed > binary_thresh, 1, 0)

            out = nib.Nifti1Image(transformed, static.affine,
                                  header=static.header)
            out.to_filename(output_filename)

    else:

        static_mask = nib.load(static_mask_path).get_fdata()

        map_mni_to_subj = find_transform(mni_fa_path, fa_path,
                                         hard_static_mask=static_mask)

        for roi_path in os.listdir(rois_path):

            if '.nii.gz' not in roi_path:
                continue

            output_filename = output_path+roi_path.split('/')[-1]

            apply_transform(rois_path+roi_path, map_mni_to_subj,
                            static_file=fa_path,
                            output_path=output_filename, binary=True,
                            binary_thresh=binary_thresh)


def tract_to_trk(input_file: str, space_file: str):
    '''


    Parameters
    ----------
    input_file : str
        DESCRIPTION.
    output_file : str
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    tract = load_tractogram(input_file, space_file)

    sft_reg = StatefulTractogram(tract.streamlines, nib.load(space_file),
                                 Space.RASMM)

    save_tractogram(sft_reg, input_file[:-3]+'trk', bbox_valid_check=False)


def tracking(patient: str, root: str, roi: str, side=None,
             number: int = 1000, angle: int = 15, num_rois: int = 3,
             max_length: int = 200, max_attempts: int = 2000,
             cutoff: float = 0.1, remove_outlier_dir: bool = False,
             decuss: bool = False):

    fod_file = root+'subjects/'+patient + \
        '/dMRI/ODF/MSMT-CSD/'+patient+'_MSMT-CSD_WM_ODF.nii.gz'
    if side == 'cc_center':
        tck_file = (root+'subjects/'+patient+'/dMRI/tractography/tois/' +
                    patient+'_'+region+'.tck')
    else:
        tck_file = (root+'subjects/'+patient+'/dMRI/tractography/tois/' +
                    patient+'_'+region+'_'+side+'.tck')
    if decuss:
        tck_file = tck_file[:-4]+'_decuss.tck'

    if not os.path.isdir(root+'/subjects/'+patient+'/dMRI/tractography/tois/'):
        os.mkdir(root+'/subjects/'+patient+'/dMRI/tractography/tois/')

    seed_side = side
    if decuss:
        target_side = {'left': 'right', 'right': 'left'}.get(side)
        exclude_file = output_roi_path+'cc_slice.nii.gz'
    else:
        target_side = side

    if side == 'cc_center':
        seed_file = output_roi_path+region+'.nii.gz'
        include1_file = output_roi_path+'cc_right.nii.gz'
        include2_file = output_roi_path+'cc_left.nii.gz'
    elif side == 'center':
        seed_file = output_roi_path+region+'_2_left.nii.gz'
        include1_file = output_roi_path+region+'_2_right.nii.gz'
    elif side is not None:
        seed_file = output_roi_path+region+'_1_'+seed_side+'.nii.gz'
        include1_file = output_roi_path+region+'_2_'+target_side+'.nii.gz'
        include2_file = output_roi_path+region+'_3_'+target_side+'.nii.gz'

    params = {'fod': 'msmt-CSD', 'algo': 'IFOD2', 'number': number,
              'angle': angle, 'num_rois': num_rois,
              'max_length': max_length, ' max_attempts':  max_attempts,
              'cutoff': cutoff, 'remove_outlier_dir': remove_outlier_dir}

    cmd = ('tckgen ' +
           fod_file + ' ' +
           tck_file + ' ' +
           '-include '+include1_file + ' ' +
           '-seed_image '+seed_file + ' ' +
           # '-max_attempts_per_seed ' + str(max_attempts) + ' ' +
           '-seeds '+str(number*max_attempts)+' ' +
           '-cutoff '+str(cutoff)+' ' +
           '-maxlength ' + str(max_length) + ' ' +
           '-select '+str(number)+' ' +
           '-angle '+str(angle)+' -force')

    if num_rois == 3:

        cmd += ' -include '+include2_file

    if decuss:

        cmd += ' -exclude '+exclude_file

    os.system(str(cmd))

    tract_to_trk(tck_file, fod_file)

    trk_file = tck_file[:-4]+'.trk'

    try:
        point_array = extract_nodes(trk_file)
        remove_outlier_streamlines(trk_file, point_array, outlier_ratio=0,
                                   remove_outlier_dir=remove_outlier_dir)
    except:
        print('No streamlines removed due to insufficient streamlines')

    with open(trk_file[:-4]+'.txt', 'w') as outfile:
        json.dump(params, outfile)


if __name__ == '__main__':

    root = sys.argv[1]
    patient = sys.argv[2]

    fa_path = (root + 'subjects/' + patient + '/dMRI/microstructure/dti/'
               + patient + '_FA.nii.gz')
    T1_path = (root+'subjects/' + patient + '/T1/' + patient +
               '_T1_brain_mask.nii.gz')
    mni_fa_path = './atlas/FSL_HCP1065_FA_1mm.nii.gz'
    mni_T1_path = './atlas/MNI152_T1_1mm_brain.nii.gz'
    rois_path = './atlas_rois/'
    static_mask_path = (root + 'subjects/' + patient + '/masks/' + patient
                        + '_brain_mask.nii.gz')

    output_roi_path = root + 'subjects/' + patient + '/dMRI/tractography/rois/'

    if not os.path.isdir(output_roi_path):
        os.mkdir(output_roi_path)
    output_roi_path = output_roi_path+patient+'_'

    register_rois_to_subj(fa_path, rois_path, mni_fa_path,
                          mni_T1_path, T1_path,
                          output_roi_path, static_mask_path)  # ,
    # T1_based_regis=True, binary_thresh=.75)  # !!!

    for region in ['cc_ant_midbody', 'cc_post_midbody', 'cc_genu', 'cc_isthmus',
                   'cc_splenium']:
        tracking(patient, root, roi=region, side='cc_center', number=2000,
                 max_length=175)

    for region in ['uf', 'cingulum']:
        tracking(patient, root, roi=region, side='left')
        tracking(patient, root, roi=region, side='right')

    for region in ['fornix']:
        tracking(patient, root, roi=region, side='left', max_attempts=2000)
    tracking(patient, root, roi=region, side='right', max_attempts=2000)

    for region in ['atr']:
        tracking(patient, root, roi=region, side='left', number=2000)
        tracking(patient, root, roi=region, side='right', number=2000)

    # for region in ['fat']:
    # tracking(patient, root, roi=region, side='left', num_rois=2)
    # tracking(patient, root, roi=region, side='right', num_rois=2)

    # for region in ['cst']:
    #     tracking(patient, root, roi=region,
    #              side='left', num_rois=2, number=2000)
    #     tracking(patient, root, roi=region,
    #              side='right', num_rois=2, number=2000)

    for region in ['ci']:
        tracking(patient, root, roi=region, remove_outlier_dir=True,
                 side='left', num_rois=2, number=4000)
        tracking(patient, root, roi=region, remove_outlier_dir=True,
                 side='right', num_rois=2, number=4000)

    # for region in ['cl']:
    #     tracking(patient, root, roi=region, side='left', max_attempts=10000,
    #              num_rois=2, number=1500)
    #     tracking(patient, root, roi=region, side='right', max_attempts=10000,
    #              num_rois=2, number=1500)
    #     tracking(patient, root, roi=region, side='left', max_attempts=25000,
    #              num_rois=2, number=1500, decuss=True)
    #     tracking(patient, root, roi=region, side='right', max_attempts=25000,
    #              num_rois=2, number=1500, decuss=True)
    #     tracking(patient, root, roi=region, side='center',
    #              num_rois=2, number=1500)

    # for region in ['slf']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=3)
    #     tracking(patient, root, roi=region, side='right', number=2000,
    #              num_rois=3)

    # for region in ['ilf']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=2)

    # for region in ['ss']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=2)
    #     tracking(patient, root, roi=region, side='right', number=2000,
    #              num_rois=2)

    # for region in ['ifof']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=2)
    #     tracking(patient, root, roi=region, side='right', number=2000,
    #              num_rois=2)

    # for region in ['af']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=2)
    #     tracking(patient, root, roi=region, side='right', number=2000,
    #              num_rois=2)

    # for region in ['vof']:
    #     tracking(patient, root, roi=region, side='left', number=2000,
    #              num_rois=2)

    # ARC -------------------------------

    # for region in ['cst']:
    # For 42 slices -----------------------------------

    # tracking(patient, root, roi=region, side='left', num_rois=2, angle=15,
    #           number=2000, max_attempts=3000, cutoff=0.08)
    # tracking(patient, root, roi=region, side='right', num_rois=2, angle=15,
    #           number=2000, max_attempts=3000, cutoff=0.08)

    # Else -----------------------------------------------
    # tracking(patient, root, roi=region, side='left', angle=15,
    #          number=2000, max_attempts=3000, cutoff=0.08)
    # tracking(patient, root, roi=region, side='right', angle=15,
    #          number=2000, max_attempts=3000, cutoff=0.08)
