import os, sys
import genet.utils
from Bio import Entrez, GenBank, SeqIO

class GetGene:
    def __init__(self, 
                 gene_name:str,
                 species:str = 'Homo sapiens',
                 search_option:str = 'AND biomol_genomic[PROP] AND RefSeqGene[Filter]',
                 ):


        print('Find %s from NCBI nucleotide database' % gene_name)
        search_string = '%s[Gene] AND %s[title] AND %s[Organism] %s' % (gene_name, gene_name, species, search_option)
        
        self.handle      = Entrez.esearch(db="nucleotide", term=search_string)
        self.gene_record = Entrez.read(self.handle)
        self.ids         = self.gene_record['IdList']
        if len(self.gene_record['IdList']) > 1:
            print('[Warnning] There are more than one ID from result. Please check your search options.')


        print('RefGenID found: ', self.ids)
        print('')

        self.fetch      = Entrez.efetch(db='nucleotide', id=self.gene_record['IdList'], rettype='gb', retmode='xlm')
        self.seq_record = SeqIO.read(self.fetch, 'genbank') 

    def is_misc_feat(self, feat): return feat.type == 'misc_feature'
    def is_source(self, feat):    return feat.type == 'source'

    def exons(self):
        def is_exon(feat):      return feat.type == 'exon'
        self.feat = self.seq_record.features
        list_exons = [f for f in filter(is_exon, self.feat)]
        return list_exons

    def transcripts(self):
        def is_mrna(feat):      return feat.type == 'mRNA'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_mrna, self.feat)]

        return list_transcripts

    def cds(self):
        def is_cds(feat):      return feat.type == 'CDS'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_cds, self.feat)]

        return list_transcripts
    

    def misc(self):
        def is_misc_feat(feat): return feat.type == 'misc_feature'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_misc_feat, self.feat)]

        return list_transcripts

    def source(self):
        def is_source(feat):    return feat.type == 'source'

        self.feat = self.seq_record.features
        list_transcripts = [f for f in filter(is_source, self.feat)]

        return list_transcripts
    

class GetClinVar:
    def __init__(self, record_id:str):


        self._record_id = record_id

        if self._record_id.startswith('VCV'):
            self.handle = Entrez.efetch(db='clinvar', id=self._record_id.split('.')[0], rettype='vcv') # VCV로 받을 경우    
        else:            
            self.handle = Entrez.efetch(db='clinvar', id=self._record_id, rettype='vcv', is_varationid='true', from_esearch="true") # variation ID로 받을 경우
        
        import xml.etree.ElementTree as ET
        self.result = ET.parse(self.handle)
        self.root = self.result.getroot()
        
        self.var_loc = self.root.findall('./VariationArchive/InterpretedRecord/SimpleAllele/Location/SequenceLocation')

        for self.info in self.var_loc:
            if self.info.attrib['Assembly'] == 'GRCh38':
                self.chr_acc = self.info.attrib['Accession']
                self.start   = int(self.info.attrib['start'])
                self.stop    = int(self.info.attrib['stop'])
                self.ref_nt  = self.info.attrib['referenceAlleleVCF']
                self.alt_nt  = self.info.attrib['alternateAlleleVCF']
                self.alt_len = int(self.info.attrib['variantLength'])
                break

        if   len(self.ref_nt) == len(self.alt_nt): self.alt_type = 'sub'
        elif len(self.ref_nt) <  len(self.alt_nt): self.alt_type = 'ins'
        elif len(self.ref_nt) >  len(self.alt_nt): self.alt_type = 'del'
    
    def seq(self, context:int = 60):
        self.chr_seq_fetch = Entrez.efetch(db="nucleotide", 
                                           id=self.chr_acc, 
                                           rettype="fasta", 
                                           strand=1, 
                                           seq_start = self.start-context, 
                                           seq_stop  = self.stop+context+self.alt_len
                                           )

        self.ref_seq = str(SeqIO.read(self.chr_seq_fetch, "fasta").seq)
        self.chr_seq_fetch.close()
        
        if self.alt_type != 'del':
            self.alt_seq = self.ref_seq[:context] + self.alt_nt + self.ref_seq[context+1:]
        else:
            self.alt_seq = self.ref_seq[:context] + self.ref_seq[context+self.alt_len:]

        if self.alt_type == 'ins':
            self.ref_seq = self.ref_seq[1:]
            self.alt_seq = self.alt_seq[1:]

        return self.ref_seq[:1+context*2], self.alt_seq[:1+context*2]