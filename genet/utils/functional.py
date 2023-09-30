import os

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower) for k, v in input.items())

class SplitFastq:
    # A function that evenly divides a fastq file into the desired number of parts.
    def __init__(
        self,
        file:str, # file (str): Path to the fastq file
        n_split:int, # n_split (int): The number of parts to divide the file into
        out_name:str, # out_name (str): Prefix for the files that will be saved after splitting
        out_path:str='./', # out_path (str, optional): The path where the output will be saved. Defaults to './'.
        silence:bool=False, # silence (bool, optional): Used to turn off print messages for logging. Defaults to False.
        ):
        
        output_format = 'fastq'
        lineset = 4

        self.names = []
        self.dir   = '%s/%s_subsplits' % (os.path.abspath(out_path), out_name)
        os.makedirs(self.dir, exist_ok=True)

        with open(file, 'r') as f:
            lines   = f.readlines()
            total   = len(lines)
            rec_cnt = total / lineset

            list_nBins = [[int(rec_cnt * (i + 0) / n_split), int(rec_cnt * (i + 1) / n_split)] for i in range(n_split)]
            self.meta  = {}
            cnt = 0

            for nStart, nEnd in list_nBins:
                if silence == False: print('[Info] Make data subsplits: %s - %s' % (nStart, nEnd))

                sSplit_file_name = '%s_%s.%s' % (out_name, cnt, output_format)
                with open('%s/%s' % (self.dir, sSplit_file_name), 'w') as outfile:
                    for l in lines[nStart*lineset:nEnd*lineset]: outfile.write(l)
                
                self.names.append(sSplit_file_name)
                
                
                self.meta[sSplit_file_name] = {
                    'start': nStart,
                    'end'  : nEnd,
                    'count': nEnd - nStart
                }
                cnt += 1