import os

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower) for k, v in input.items())

class SplitFastq:
    def __init__(
        self,
        file:str,
        n_split:int,
        out_name:str,
        out_path:str='./',
        silence:bool=False,
        ):
        """fastq file을 원하는 수 만큼 균등하게 나눠주는 함수.

        Args:
            file (str): fastq 파일 경로
            n_split (int): 몇 등분 할 것인지 적는 칸
            out_name (str): 나눴을 때 저장되는 파일들의 prefix
            out_path (str, optional): Output이 저장 될 경로. Defaults to './'.
            silence (bool, optional): Logging을 위한 print 되는 메시지를 끄는 용도. Defaults to False.
        """        
        
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


# class END: SplitFastq