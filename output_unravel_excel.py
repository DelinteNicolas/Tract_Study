# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:27:03 2024

@author: N. Delinte
"""

import sys
import json
import xlsxwriter


def dicToExcel(totDic: dict, Excel_path: str, st_luc_out: bool = False,
               verbose: bool = True, unwanted_list: list = [],
               unwanted_metrics: list = [], unwanted_regions: list = []):

    workbook = xlsxwriter.Workbook(
        Excel_path+'.xlsx', {'nan_inf_to_errors': True})
    worksheet = workbook.add_worksheet('Means and Standard deviations')

    Mean_dic = totDic['Mean']
    try:
        Dev_dic = totDic['Dev']
        dev = True
    except KeyError:
        dev = False

    PatientList = list(Mean_dic.keys())
    PatientList.sort()
    for p in unwanted_list:
        try:
            PatientList.remove(p)
        except ValueError:
            continue
    RegionList = list(set().union(*Mean_dic.values()))
    RegionList.sort()
    for r in unwanted_regions:
        try:
            RegionList.remove(r)
        except ValueError:
            continue
    MetricList = set()
    for p in PatientList:
        MetricList.update(set().union(*Mean_dic[p].values()))
    MetricList = list(MetricList)
    MetricList.sort()
    MetricList = list(filter(lambda x: 'DIFF' not in x, MetricList))
    for m in unwanted_metrics:
        try:
            MetricList.remove(m)
        except ValueError:
            continue

    for p in PatientList:
        for r in RegionList:
            try:
                Mean_dic[p][r][MetricList[0]]
            except KeyError:
                print('No region called '+r+' in patient '+p)
                continue

            if st_luc_out:
                dec = 1
            else:
                dec = 0

            for i in MetricList:

                try:

                    worksheet.write(MetricList.index(i)+(len(MetricList)+2)
                                    * PatientList.index(p)+1,
                                    RegionList.index(r)+1+dec,
                                    Mean_dic[p][r][i])
                    if dev:
                        worksheet.write(MetricList.index(i)+(len(MetricList)+2)
                                        * PatientList.index(p)+1,
                                        RegionList.index(r) +
                                        1+dec+len(RegionList)+2+dec+dec,
                                        Dev_dic[p][r][i])

                except KeyError:
                    if verbose:
                        print(i+' not found for '+p+' : '+r)

            if st_luc_out:

                worksheet.write((len(MetricList)+2)*PatientList.index(p),
                                RegionList.index(r)+1+dec,
                                r.replace('.', '_'))
                if dev:
                    worksheet.write((len(MetricList)+2)*PatientList.index(p),
                                    RegionList.index(r)+1+dec +
                                    len(RegionList)+dec
                                    + 2+dec, r.replace('.', '_'))

            else:

                worksheet.write((len(MetricList)+2)*PatientList.index(p),
                                RegionList.index(r)+1+dec,
                                p+'_'+r.replace('.', '_'))
                if dev:
                    worksheet.write((len(MetricList)+2)*PatientList.index(p),
                                    RegionList.index(r)+1+dec+len(RegionList)
                                    + 2+dec, p+'_'+r.replace('.', '_'))

        if st_luc_out:

            for i in MetricList:
                worksheet.write(list(MetricList).index(i)+(len(MetricList)+2)
                                * PatientList.index(p)+1, 0, i)
                worksheet.write(list(MetricList).index(i)+(len(MetricList)+2)
                                * PatientList.index(p)+1, 1, p)
                if dev:
                    worksheet.write(list(MetricList).index(i)+(len(MetricList)+2)
                                    * PatientList.index(p)+1,
                                    len(RegionList)+4, i)
                    worksheet.write(list(MetricList).index(i)+(len(MetricList)+2)
                                    * PatientList.index(p)+1,
                                    len(RegionList)+5, p)

        else:
            for i in MetricList:
                worksheet.write(MetricList.index(
                    i)+(len(MetricList)+2)*PatientList.index(p)+1, 0, i)
                worksheet.write(MetricList.index(i)+(len(MetricList)+2)
                                * PatientList.index(p)+1,
                                len(RegionList)+2, i)

    workbook.close()


if __name__ == '__main__':

    st_luc_out = True
    dic_file = sys.argv[1]+'out/unravel_mean_ang_tsl.json'
    unwanted_metrics = ['snr', 'fintra', 'fextra', 'fiso', 'odi']
    unwanted_regions = []
    unwanted_list = []

    excel_file = dic_file[:-5]

    # file = open(unwanted_list_path)
    # unwanted_list = json.load(file)
    # file.close()

    tot = json.load(open(dic_file,))
    dicToExcel(tot, excel_file, st_luc_out=st_luc_out,
               unwanted_list=unwanted_list, unwanted_metrics=unwanted_metrics,
               unwanted_regions=unwanted_regions)
