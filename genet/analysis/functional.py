import genet
import os, sys, regex, glob, shutil, itertools, time, subprocess
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO
from tqdm import tqdm

class SortByBarcodes:
    def __init__(self,
                 seq_file:str,
                 barcode_file:str,
                 barcode_pattern:str = None,
                 output_name:str = 'barcode_sorted', 
                 output_path:str = './',
                 data_format:str = 'fastq',
                 output_format:str = 'fastq',
                 n_cores:int = int(mp.cpu_count()*0.5),
                 remove_temp_files:bool = True,
                 silence:bool = False,
                 ):

        # check input types
        if n_cores > mp.cpu_count():
            sys.exit('n_core should be lower than the number of cores which your machine has')

        # load barcode and data files
        self.df_bc    = pd.read_csv(barcode_file, names=['id', 'barcode'])

        splits = genet.utils.SplitFastq(
            seq_file, n_cores, out_path=output_path, out_name=output_name, silence=silence)
        
        self.list_sParameters = []
        for s in splits.names:
            self.list_sParameters.append([self.df_bc, barcode_pattern, s, splits.dir, output_format, silence])

        # sorting barcodes from subsplit files
        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: sorting by barcodes')
        p.map_async(sort_barcode, self.list_sParameters).get()

        p.close()
        p.join()

        # combine all temp files
        if silence == False: print('[Info] Make final sorted files')
        sOUT_DIR = '%s/%s_sorted' % (output_path, output_name)
        os.makedirs(sOUT_DIR, exist_ok=True)

        self.couts = mp.Manager().dict()
        self.couts['_Not_matched'] = []
        for key in self.df_bc['barcode']: self.couts[key] = []

        self.barcodes = self.couts.keys()
        self.nBarcodeCnt = len(self.barcodes)

        list_barcode_nBins = [[int(self.nBarcodeCnt * (i + 0) / n_cores), int(self.nBarcodeCnt * (i + 1) / n_cores)] for i in range(n_cores)]

        if silence == False: print('[Info] Make barcode subsplits')
        self.list_combine_param = [[splits.dir, output_format, sOUT_DIR, self.couts, self.barcodes[nStart:nEnd], silence] for nStart, nEnd in list_barcode_nBins]

        p = mp.Pool(n_cores)
        if silence == False: print('[Info] Starting map_async: combine all temp files')
        p.map_async(combine_files, self.list_combine_param).get()

        p.close()
        p.join()

        # finalize
        if remove_temp_files==True:
            if silence == False: print('[Info] Removing temp files')
            shutil.rmtree(splits.dir)
        if silence==False: print('[Info] Done: SortByBarcodes - %s' % output_name)


def sort_barcode(list_sParameters):
    df_bc            = list_sParameters[0]
    barcode_pattern  = list_sParameters[1]
    subsplit_name    = list_sParameters[2]
    subsplit_dir     = list_sParameters[3]
    output_format    = list_sParameters[4]
    silence          = list_sParameters[5]
    
    if silence == False: print('[Info] Barcode sorting - %s' % (subsplit_name))
    
    # make temp dir
    temp_dir = '%s/_temp_%s_sorting' % (subsplit_dir, subsplit_name)
    os.makedirs(temp_dir, exist_ok=True)

    dict_barcode = {'_Not_matched': []}
    for key in df_bc['barcode']:
        dict_barcode[key] = []
    
    fq_file = '%s/%s' % (subsplit_dir, subsplit_name)
    record_iter = SeqIO.parse(open(fq_file), output_format)

    for rec in record_iter:
        seq = str(rec.seq)
        check_match = False
        
        if barcode_pattern == None:
            for k in dict_barcode.keys():
                if k not in seq: continue
                else:
                    dict_barcode[k].append(rec)
                    check_match=True
                    break

        else:
            try:
                for sReIndex in regex.finditer(barcode_pattern, seq, overlapped=True):
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end()
                    window = seq[nIndexStart:nIndexEnd]
                    
                    try: 
                        dict_barcode[window].append(rec)
                        check_match = True
                        break
                    except KeyError: continue
            except KeyError: continue
        
        if check_match==False: dict_barcode['_Not_matched'].append(rec)
            
    if silence == False: print('Make temp sorted %s file: %s' % (output_format, subsplit_name))
    
    for barcode, seq_rec in dict_barcode.items():
        SeqIO.write(seq_rec, '%s/%s.%s' % (temp_dir, barcode, output_format), output_format)


def combine_files(list_combine_param):
    splits_dir    = list_combine_param[0]
    output_format = list_combine_param[1]
    sOUT_DIR      = list_combine_param[2]
    counts        = list_combine_param[3]
    barcodes      = list_combine_param[4]
    silence       = list_combine_param[5]

    for key in barcodes:
        if silence == False: print('Make combined file: %s' % key)
        temp_fqs = glob.glob('%s/**/%s.%s' % (splits_dir, key, output_format))

        output_file_name = '%s/%s.%s' % (sOUT_DIR, key, output_format)

        with open(output_file_name, 'w') as outfile:
            for filename in sorted(temp_fqs):
                with open(filename) as file:        
                    outfile.write(file.read())
        
        with open(output_file_name, 'r') as outfile:
            counts[key] = len(outfile.readlines())


def sort_barcode_and_combine(list_sParameters):
    df_bc            = list_sParameters[0]
    barcode_pattern  = list_sParameters[1]
    subsplit_name    = list_sParameters[2]
    subsplit_dir     = list_sParameters[3]
    output_format    = list_sParameters[4]
    silence          = list_sParameters[5]
    
    if silence == False: print('[Info] Barcode sorting - %s' % (subsplit_name))
    
    # make temp dir
    temp_dir = '%s/_temp_%s_sorting' % (subsplit_dir, subsplit_name)
    os.makedirs(temp_dir, exist_ok=True)

    dict_barcode = {'_Not_matched': []}
    for key in df_bc['barcode']:
        dict_barcode[key] = []
    
    fq_file = '%s/%s' % (subsplit_dir, subsplit_name)
    record_iter = SeqIO.parse(open(fq_file), output_format)

    for rec in record_iter:
        seq = str(rec.seq)
        check_match = False
        
        if barcode_pattern == None:
            for k in dict_barcode.keys():
                if k not in seq: continue
                else:
                    dict_barcode[k].append(rec)
                    check_match=True
                    break

        else:
            try:
                for sReIndex in regex.finditer(barcode_pattern, seq, overlapped=True):
                    nIndexStart = sReIndex.start()
                    nIndexEnd = sReIndex.end()
                    window = seq[nIndexStart:nIndexEnd]
                    
                    try: 
                        dict_barcode[window].append(rec)
                        check_match = True
                        break
                    except KeyError: continue
            except KeyError: continue
        
        if check_match==False: dict_barcode['_Not_matched'].append(rec)
            
    if silence == False: print('Make temp sorted %s file: %s' % (output_format, subsplit_name))
    
    for barcode, seq_rec in dict_barcode.items():
        SeqIO.write(seq_rec, '%s/%s.%s' % (temp_dir, barcode, output_format), output_format)