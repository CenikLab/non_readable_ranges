from ribopy import Ribo
import ribopy
from Fasta import FastaFile
import json
from argparse import ArgumentParser

study = 'GSE51584'
experiment = 'GSM1248729'
default_ribo_path =  f"/scratch/users/mjgeng/process-multiple-ribo/output/{study}/ribo/experiments/{experiment}.ribo"

parser = ArgumentParser()
parser.add_argument('--path_to_fasta', default="data/appris_human_v2_selected.fa.gz", dest='path_to_fasta')
parser.add_argument('--path_to_ribo', default=default_ribo_path, dest='path_to_ribo')
args = parser.parse_args()

ribo = Ribo(args.path_to_ribo, alias=ribopy.api.alias.apris_human_alias)

fasta = FastaFile(args.path_to_fasta)
fasta_dict = {e.header: e.sequence for e in fasta}
sequence_dict = {
    ribopy.api.alias.apris_human_alias(transcript): fasta_dict[transcript] for transcript in ribo.transcript_names
}

with open("data/sequence_dict.json", 'w+') as f:
    json.dump(sequence_dict, f)