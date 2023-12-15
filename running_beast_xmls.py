# -*- coding: utf-8 -*-
# @Author: Martin Grunnill
# @Date:   2023-12-04 20:09:41
# @Last Modified by:   Martin Grunnill
# @Last Modified time: 2023-12-06 10:36:42
import os 
import shutil

working_dir = os.getcwd()
beast_xmls = os.path.join(working_dir, 'beast_xmls')
data_sets = [#'ON-2022','MD-2022',
             'ON-under5s-2022',
             'n20-ON-under5s-2022', 'n30-ON-under5s-2022',
             'n40-ON-under5s-2022','n50-ON-under5s-2022'             
             ]

for data_set in data_sets:
    data_set_dir = os.path.join(beast_xmls, data_set)
    data_set_xml = data_set + '.xml'
    for run in range(7,17):
        run_name = 'run_'+str(run)
        run_dir = os.path.join(data_set_dir, run_name)
        run_xml = run_name+'.xml'
        os.mkdir(run_dir)        
        shutil.copyfile(os.path.join(data_set_dir,data_set_xml), os.path.join(run_dir,run_xml))
        os.chdir(run_dir)
        os.system('gnome-terminal -- beast -beagle_SSE ' + run_xml)
        os.chdir(working_dir)
        