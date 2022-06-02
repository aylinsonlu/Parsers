import subprocess
import os
from pprint import pprint
import sys

def get_cds(gtf_file,gene,genome_file):
    cds_dict = {}
    with open(gtf_file,'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.rstrip().split('\t')
                attributes = fields[8].replace('"',"")
                gene_id = attributes.split("gene_id ")[1].split(";")[0]
                type = fields[2]
                chrm = fields[0]
                strand = fields[6]
                if gene_id == gene:
                    if type == "CDS":
                        exon =  attributes.split("exon_number ")[1].split(";")[0]
                        start = fields[3]
                        stop = fields[4]
                        transcript_id = attributes.split("transcript_id ")[1].split(";")[0]
                        if strand=="-":                            
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop + " -i").read().split("\n")
                            seq = "".join(output[1:])
                        if strand=="+":
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop).read().split("\n")
                            seq = "".join(output[1:])
                        if not cds_dict:
                            cds_dict[gene_id] = {'chrm':chrm,'strand':strand,'transcripts':{transcript_id:[{'exon':exon,'start':start,'stop':stop,'seq':seq}]}}
                        else:
                            if transcript_id not in cds_dict[gene_id]['transcripts'].keys():
                                cds_dict[gene_id]['transcripts'].update({transcript_id:[{'exon':exon,'start':start,'stop':stop,'seq':seq}]})
                            else:
                                cds_dict[gene_id]['transcripts'][transcript_id].append({'exon':exon,'start':start,'stop':stop,'seq':seq})
                        
    print(cds_dict)
    return cds_dict


