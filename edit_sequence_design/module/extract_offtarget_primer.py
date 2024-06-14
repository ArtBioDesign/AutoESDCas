import os
import pandas as pd
import numpy as np 
import warnings          
warnings.filterwarnings('ignore')  
import re  
from os.path import exists
import sgRNA_utils.sgRNA_primer_util as util  
from sgRNA_utils.sgRNA_primer_config import config  


def get_invalid_primer(invalid_priemr_df):
    li =[]
    for i,v in invalid_priemr_df.iterrows():  
            
            state,id,cor,name,offtarget = v["ID"].split(":")  #PCR
            id = id + ":" + cor                               #PCR
            off_target = str(v["off_target"])
            off_target = eval(off_target)

            print(id)   #Pseudomonas_aeruginosa 500-1
    
            #目标引物  
            SEQUENCING_PRIMER_1 = ""
            SEQUENCING_PRIMER_1_offtarget = ""
            SEQUENCING_PRIMER_2 = ""
            SEQUENCING_PRIMER_2_offtarget = ""

            #脱靶引物
            coordinates_1 = []  
            coordinates_2 = []
            for i,v in off_target.items():
                
                if len(v.split(";")) == 1 and "SEQUENCING_PRIMER_1" in i:
                    SEQUENCING_PRIMER_1_0 = v    
                    chrom_cor_strand = SEQUENCING_PRIMER_1_0.split(";")[0]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_1 =  chrom_cor
                
                if  len(v.split(";")) == 1 and "SEQUENCING_PRIMER_2" in i:
                    SEQUENCING_PRIMER_2_0 = v
                    chrom_cor_strand = SEQUENCING_PRIMER_2_0.split(";")[0]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_2 = chrom_cor 
            
            # for i,v in off_target.items(): 
                chrom_cor_strand = v.split(";")[0]
                chrom_cor,strand = chrom_cor_strand.split("|")

                if "SEQUENCING_PRIMER_1" in i and len(v.split(";")) != 1:  
                    coordinates_1.append( chrom_cor )
                if "SEQUENCING_PRIMER_2" in i and len(v.split(";")) != 1: 
                    coordinates_2.append( chrom_cor )
            
            off_target_type = "" 
            if  SEQUENCING_PRIMER_1 !="": 
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)  
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" 
            elif  SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_2 offtarget"
            if SEQUENCING_PRIMER_1 !="" and SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" + "&" +  "SEQUENCING_PRIMER_2 offtarget"
                
            dict_offtarget = {
                    "locus_tag" : id.split(";")[0],  
                    "SEQUENCING_PRIMER_1 location (start-end)":SEQUENCING_PRIMER_1,
                    "SEQUENCING_PRIMER_1 off-target site":SEQUENCING_PRIMER_1_offtarget,
                    "SEQUENCING_PRIMER_2 location (start-end)":SEQUENCING_PRIMER_2,
                    "SEQUENCING_PRIMER_2 off-target site":SEQUENCING_PRIMER_2_offtarget,
                    "type (B: Before optimization);(A:After optimization)":off_target_type
                }
            li.append(dict_offtarget)   
    
    if len(li)>0:
        df = pd.DataFrame(li)
    else:
        df = pd.DataFrame(columns=["locus_tag"])
    
    return df

def extract_failture_primer(one_plasmid_design_F_Result_id):

    fail_primer_df = pd.DataFrame( columns=["locus_tag"] ,data = list(one_plasmid_design_F_Result_id) )  
    fail_primer_df["locus_tag"] = fail_primer_df.locus_tag.apply(lambda x: x.split(";")[0] )
    fail_primer_df["type (B: Before optimization);(A:After optimization)"] = "No primer3 design results"
    
    return fail_primer_df

def df_to_fasta(df,fasta_filename,id,sequence_name):
    with open(fasta_filename, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence_id = row[id]
            sequence = row[sequence_name]
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')

def rename_invalid_primer(invalid_primer, opt_type ):

    invalid_primer = invalid_primer.rename( columns={
                                                        "SEQUENCING_PRIMER_1 location (start-end)": f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)",
                                                        "SEQUENCING_PRIMER_1 off-target site":f"{opt_type} optimization SEQUENCING_PRIMER_1 off-target site",
                                                        "SEQUENCING_PRIMER_2 location (start-end)":f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)",
                                                        "SEQUENCING_PRIMER_2 off-target site":f"{opt_type} optimization  SEQUENCING_PRIMER_2 off-target site"
                                                    } )
    
    if  opt_type == "Before":
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].apply(lambda x: "B:"+x)
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].astype(str)
    elif opt_type == "After":
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].astype(str)
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].apply(lambda x:  "All primer templates offtarget"  if x == "" else x )
        invalid_primer["type (B: Before optimization);(A:After optimization)"] =  "A:" + invalid_primer["type (B: Before optimization);(A:After optimization)"]

    
    return invalid_primer

def extract_opt_invalid_primer(Test_primer_G, invalid_primer, genome_path, opt_type, opt_type_value = "A:All primer templates offtarget"):

    ab = Test_primer_G.copy()  
    ab["ID"] = list( ab.ID.apply(lambda x: x.split(";")[0]) )
    ab.index = ab["ID"]

    for i,v in invalid_primer.iterrows():
        id = v["locus_tag"]
        if v["type (B: Before optimization);(A:After optimization)"] != opt_type_value:
                
                if pd.isna( v[f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] ) or  v[f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] == "":
                    left_primer = ab.loc[id,"SEQUENCING_PRIMER_1"]
                    seq_df = pd.DataFrame(
                            [{
                            "id":id,
                            "sequence":left_primer
                            }]
                        )
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_1", genome_path, "0")
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_1 location (start-end)"]

                if pd.isna(  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]  ) or  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]=="":
                    

                    right_primer = ab.loc[id,"SEQUENCING_PRIMER_2"]
                    
                    seq_df = pd.DataFrame(
                            [{
                            "id":id,
                            "sequence":right_primer
                            }]
                            
                        )
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_2", genome_path, "16")
            
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_2 location (start-end)"]

    return invalid_primer

def parse_sam_file(sam_file_path = "alignment.sam"):

    # 定义用于存储比对结果的列表
    alignment_data = []  

    with open(sam_file_path, "r") as sam_file:
        for line in sam_file:
            if not line.startswith("@"):  # 跳过SAM文件头部
                fields = line.strip().split("\t")
                read_name = fields[0]  
                chain = fields[1]
                reference_name = fields[2]
                reference_start = int(fields[3])  
                sequence = fields[9]
                mismatch = fields[-2]
                matching_number = fields[-1]

                alignment_data.append([read_name, chain, reference_name, reference_start, sequence, mismatch, matching_number])

    # 创建DataFrame
    columns = ["ReadName","chain", "ReferenceName", "ReferenceStart", "Sequence","Mismatch", "MatchingNumber"]
    alignment_df = pd.DataFrame(alignment_data, columns=columns)

    return alignment_df

def get_no_offtarget(no_off_target,name):

    if len( no_off_target ) > 0:
        no_off_target["len"] = no_off_target.Sequence.apply(lambda x: len(x))
        no_off_target["ReferenceStart"] = no_off_target["ReferenceStart"]
        no_off_target["ReferenceEnd"] = no_off_target["ReferenceStart"] + no_off_target["len"]
        no_off_target["ReferenceStart"]  = no_off_target["ReferenceStart"].astype("str")
        no_off_target["ReferenceEnd"] = no_off_target["ReferenceEnd"].astype("str")
        no_off_target[f"After optimization {name} location (start-end)"] = no_off_target["ReferenceName"] +":"+ no_off_target["ReferenceStart"]+"-"+no_off_target["ReferenceEnd"]
    else:
        no_off_target = pd.DataFrame( columns = ["ReadName", f"After optimization {name} location (start-end)"]  )

    return no_off_target

def algin(A_no_offtarget_B_no_result_primer_df,sequence_name, genome_path, chain):

    index_prefix_path = os.path.join(config.tmp_path, 'genome_index')
    fasta_filename = os.path.join( config.tmp_path, 'sequences.fasta')  

    util.df_to_fasta(A_no_offtarget_B_no_result_primer_df, fasta_filename)  
    sam_df = util.bowtie_seq_genome(config.tmp_path, genome_path, index_prefix_path, fasta_filename)

    # no_off_target = sam_df[(sam_df['MatchingNumber'] == 'XM:i:1') & (sam_df['Mismatch'] == 'NM:i:0') ] 
    temp = sam_df.groupby("ReadName").apply(lambda x: x[x["chain"]==chain] if len( x[x["chain"]==chain]) == 1 else None )    

    no_off_target = temp.reset_index(drop=True)
    no_off_target = get_no_offtarget(no_off_target, sequence_name)
    no_off_target = no_off_target[["ReadName", f"After optimization {sequence_name} location (start-end)"]]
    no_off_target = no_off_target.rename(columns={"ReadName":"ID"})

    return no_off_target


