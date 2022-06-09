import subprocess
import os
from pprint import pprint
import sys

def get_all_isoform_cds(gtf_file,genome_file):
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
                if type == "CDS":
                    exon =  attributes.split("exon_number ")[1].split(";")[0]
                    start = fields[3]
                    stop = fields[4]
                    transcript_id = attributes.split("transcript_id ")[1].split(";")[0]
                    gene_bio_type = attributes.split("gene_biotype ")[1].split(";")[0]
                    if gene_bio_type == "protein_coding":
                        if strand=="-":                            
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop + " -i").read().split("\n")
                            seq = "".join(output[1:])
                        if strand=="+":
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop).read().split("\n")
                            seq = "".join(output[1:])
                        if not cds_dict:
                            cds_dict[gene_id] = {'chrm':chrm,'strand':strand,'transcripts':{transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]}}
                        else:
                            if gene_id in cds_dict.keys():
                                if transcript_id not in cds_dict[gene_id]['transcripts'].keys():
                                    cds_dict[gene_id]['transcripts'].update({transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]})
                                else:
                                    cds_dict[gene_id]['transcripts'][transcript_id].append({'exon':exon,'coordinates':(start,stop),'seq':seq})
                            else:
                                cds_dict[gene_id] = {'chrm':chrm,'strand':strand,'transcripts':{transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]}}

    #print(cds_dict)
    return cds_dict

def get_single_transcript_cds(cds_dict,transcript_id):
    cds = ""
    for gene in cds_dict.keys():
        for transcript in cds_dict[gene]['transcripts'].keys():
            if transcript == transcript_id:
                for exons in cds_dict[gene]['transcripts'][transcript_id]:
                    cds +=  exons['seq']
                  
    
  
    return cds.upper()


