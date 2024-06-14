#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :test_time.py
# @Time         :2024/01/27 14:38:12   
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content

import random  
import pandas as pd
from Bio import SeqIO  
from os.path import exists,dirname,splitext,basename,realpath,abspath
import os,sys
from pprintpp import pprint as pp  
import multiprocessing  
from Bio import SeqIO  
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
workdir_tmp="/home/yanghe/tmp/"

def blast_search(input_file_path,genome,workdir):
    """
    Arguments:
        input_file_path[str]:  input information
        genome[str]:  edit reference
        workdir[str]: dir path
    Return[dict]
        {
            "key1":{
                "chrom":"NC_001133.9",
                "start":"1000",
                "unique_mapped":1,  100% comparison
                "mapped":1,  frequency
                "reverse_mapped":1, frequency
                "description":"NC_001133.9:26529-26708;NC_001133.9:26664-26843;",  
            },
        }
    """
    blast_output_file_path=os.path.join(workdir,'blast_search.txt')
    ref_lib=genome.split('/')[-1].split('.')[0]
    ref_lib = os.path.join(workdir, ref_lib)

    input_fasta = os.path.join(workdir,'blast.fasta')
    fasta_length_dict = {}
    my_records = []
    seqid=0
    with open(input_file_path,'r') as ifile:
        index = 0 
        for line in ifile:
            if index != 0 :
                seq = line
                seq_len= len(seq)-1
            index += 1
    # input fasta
    SeqIO.write(my_records, input_fasta, "fasta")
    
    # run blast
    cmd1 = "makeblastdb -in "+genome+" -dbtype nucl -parse_seqids -out "+ref_lib
    print(cmd1)
    os.system(cmd1)
    cmd2 = "blastn -query "+input_file_path+" -db "+ref_lib+" -outfmt 6 -task blastn -out "+blast_output_file_path+" -evalue 1e-30"
    os.system(cmd2)
    os.system("rm %s.n*" % ref_lib)  

    # return
    dictall = -1
    with open(blast_output_file_path,"r") as f:
        for i in f:
            linelist = i.split('\t')
            key = linelist[0]
            chrom = linelist[1]
            identity = linelist[2]
            allength = linelist[3]
            start = linelist[8]
            end = linelist[9]

            if (int(float(identity)) == 100 and int(float(allength))==seq_len):
                dictall = 1 #成功
                print("成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功成功")
            else:
                dictall = 0 ##失败
                print("失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败失败")

    return dictall




def inspect_off_target(sequence_id, seq, genome_path):

    df = pd.DataFrame(
                [{
                    "id": sequence_id,
                    "sequence": seq
                }]
            )
    input_file_path = os.path.join(workdir_tmp, "seq.fasta") 
    df_to_fasta(df, input_file_path)
    detial = blast_search(input_file_path, genome_path, workdir_tmp)
    return detial



def df_to_fasta(df,fasta_filename):
    with open(fasta_filename, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence_id = row['id']
            sequence = row['sequence']
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')


def generate_random_SNP(genome, num_snps):
    """
    生成随机的单核苷酸变异（SNP）
    参数：
    genome: 大肠杆菌基因组序列（字符串）
    num_snps: 要生成的SNP数量
    返回值：
    变异后的基因组序列（字符串）
    """
    genome_list = list(genome)
    genome_len = len(genome)

    mute_list=[]
    
    for _ in range(num_snps):

        position = random.randint(0, genome_len-1)  # 随机选择变异位置
        old_base = genome[position:position+1]
        new_base = random.choice(['A', 'T', 'C', 'G'])  # 随机选择新的碱基
        genome_list[position] = new_base  # 替换基因组中的碱基
        gene_mute = { 
                    'Manipulation type': 'substitution',
                    'position':position,
                    'Reference sequence':old_base,
                    'Inserted sequence':new_base  
                     }
        mute_list.append( gene_mute )

    mutated_genome = ''.join(genome_list)

    return mutated_genome, mute_list

def generate_random_indel(genome, num_indels):
    """
    生成随机的插入/缺失变异
    参数：
    genome: 大肠杆菌基因组序列（字符串）
    num_indels: 要生成的插入/缺失数量
    返回值：
    变异后的基因组序列（字符串）
    """
    genome_list = list(genome)
    genome_len = len(genome)

    mute_list=[]
    
    for _ in range(num_indels):
        position = random.randint(0, genome_len-1)  # 随机选择变异位置
        mutation_type = random.choice(['insertion', 'deletion'])  # 随机选择插入或缺失
        old_base = genome[position:position+1]

        if mutation_type == 'insertion':
            new_base = random.choice(['A', 'T', 'C', 'G'])  # 随机选择插入的新碱基
            genome_list.insert(position, new_base)  # 在指定位置插入新碱基
        else:
            new_base = '-'
            genome_list.pop(position)  # 在指定位置删除碱基

        gene_mute = { 
                    'Manipulation type': mutation_type,
                    'position':position,
                    'Reference sequence':old_base,
                    'Inserted sequence':new_base
                     }
        mute_list.append( gene_mute )

    mutated_genome = ''.join(genome_list)

    return mutated_genome, mute_list

def generate_mutated_genome(genome, sequence_id, mutation_rate, genome_path, li_dict):
    mutated_genome = list(genome)
  
    mute_list=[]
    mutation_index = random.randint(10000, len(mutated_genome) - 10000)    #随机生成位置

    for _ in range(int(mutation_rate )):
        mutation_type = random.choice(['deletion', 'insertion', 'substitution'])
       
        j = 0

        while True:
            # blasta序列检查
            Sequence_upstream_of_the_manipulation_site_seq = genome[mutation_index-7000 : mutation_index]

            dataill = inspect_off_target(sequence_id, Sequence_upstream_of_the_manipulation_site_seq, genome_path)
            if dataill == 1:
                print("此段不脱靶", j, file=sys.stdout)
                mutation_index  = mutation_index - 1000
                break
            elif j == 80:
                error_message = "几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复几乎整个染色体和基因组重复"
                raise ValueError(error_message)
            else: 
                j = j + 1
                print("此段脱靶此段脱靶此段脱靶此段脱靶此段脱靶此段脱靶此段脱靶此段脱靶此段脱靶此段",sequence_id,mutation_index, j, file=sys.stdout)
                sequence_id = random.choices( list(li_dict.keys()) )[0]               #随机选择染色体
                genome = li_dict[sequence_id]
                mutated_genome = list(genome)
                mutation_index = random.randint(10000, len(mutated_genome) - 10000)    #随机生成位置    

        if mutation_type == 'deletion':
            deletion_length = random.randint(1, 100)  # 调整删除片段的长度范围
            new_base = '-'
            old_base = genome[mutation_index:mutation_index + deletion_length]
            position = mutation_index, mutation_index + deletion_length

            del mutated_genome[mutation_index:mutation_index + deletion_length]

        elif mutation_type == 'insertion':
            insertion_length = random.randint(1, 100)  # 调整插入片段的长度范围
            insertion = generate_random_sequence(insertion_length)
            new_base = insertion
            old_base = '-'
            mutated_genome[mutation_index:mutation_index] = insertion

            position = mutation_index, mutation_index+1


        elif mutation_type == 'substitution':
            replacement_length = random.randint(1, 100)  # 调整替换片段的长度范围
            replacement = generate_random_sequence(replacement_length)
            new_base = replacement
            old_base = genome[mutation_index: mutation_index + replacement_length]
            mutated_genome[mutation_index:mutation_index + replacement_length] = replacement
            position = mutation_index, mutation_index + replacement_length


        gene_mute = {
                    'Name': f'{sequence_id}_{_}_{mutation_type}',
                    'Sequence upstream of the manipulation site (>100bp)': genome[mutation_index-500:mutation_index],
                    'Reference sequence':old_base,
                    'Inserted sequence':new_base,
                    'Manipulation type': mutation_type,
                    'position': position
                     }
        mute_list.append(gene_mute)

        
    return ''.join(mutated_genome), mute_list

def generate_random_sequence(length):
    bases = ['A', 'T', 'C', 'G']
    sequence = [random.choice(bases) for _ in range(length)]
    return ''.join(sequence)

def generate_one_genome_mutation(input_genome_path, input_mutation_dir,num):

    li_dict = {}
    #取基因组序列
    for record in SeqIO.parse(input_genome_path, "fasta"):
        sequence_id = record.id
        sequence = record.seq
        genome = str( record.seq )
        if len(genome) > 30000:
            li_dict.update(
                {sequence_id:genome}
            )

    path = []

    #随机生成
    li =[] 
    for i in range(0,500,50): 
        mutation_rate = i + 50  # 突变次数
        for i in range(1,num):

            input_mutation_child_dir = os.path.join( input_mutation_dir,  str(mutation_rate)+ "-" +str(i) )
                                        
           
            if not exists(input_mutation_child_dir):    
                os.makedirs(input_mutation_child_dir)
                #随机取染色体
                sequence_id = random.choices( list(li_dict.keys()) )[0]
                mutated_genome, info = generate_mutated_genome(li_dict[sequence_id], sequence_id, mutation_rate, input_genome_path, li_dict)
                df = pd.DataFrame(info)
                df.to_csv(os.path.join(input_mutation_child_dir, str(mutation_rate)+ "-" +str(i)+".csv" ), index=False)
                print( os.path.join(input_mutation_child_dir, str(mutation_rate)+ "-" +str(i)+".csv")  + "创建完毕" )
            else:
                print( os.path.join(input_mutation_child_dir, str(mutation_rate)+ "-" +str(i)+".csv")  + "已经创建完毕" )
            
            path.append( input_mutation_child_dir  )

    return path
            
def generate_batch_genome_mutation(genome_name_list, input_dir):

    genome_dict = {}

    for genome_name in genome_name_list: 

        input_genome_path  = os.path.join(input_dir, genome_name)

        name,suffix = splitext(genome_name)
        input_mutation_dir =  os.path.join(input_dir, name)

        if not exists(input_mutation_dir):
            os.makedirs(input_mutation_dir)

        path = generate_one_genome_mutation(input_genome_path, input_mutation_dir,6)

        genome_dict.update({name:path})
        
    return genome_dict

def excute_data_preproccessing(input_csv,ref_genome,output_csv,scene):
    import sys
    sys.path.append("/home/yanghe/program/data_preprocessing/")
    import parse_input_to_df as data_preprocessing  
    data = {  
                "input_file_path":input_csv,
                "ref_genome":ref_genome,
                "data_preprocessing_workdir":output_csv,
                "scene":scene
            }
    a = data_preprocessing.main(data)
    print(a)

def excute_chopchop(input_csv, ref_genome, output_csv, chopchop_config="" ):
    sys.path.append("/home/yanghe/program/chopchop/chopchop/")
    import chop_main 
    data = {
                "input_file_path": input_csv,
                "ref_genome": ref_genome,
                "chopchop_workdir": output_csv, 
                "chopchop_config":{
                    "PAM": "NGG", 
                    "guideSize": 20,
                    "maxMismatches": 3,
                    "scoringMethod": "DOENCH_2014"
                }
            }   
    a = chop_main.main(data)
    print(a)

def excute_edit_sequence_design( editor_info_input_csv, sgRNA_csv, output_csv, ref_genome, gb, scene, plasmid_metod, analysis):
    import sys
    sys.path.append('/home/yanghe/program/edit_sequence_design/')
    import edit_main

    data = {     
                "analysis":analysis,  
                "chopchop_input": editor_info_input_csv,   
                "sgRNA_result_path": sgRNA_csv, 
                "edit_sequence_design_workdir": output_csv,  
                "ref_genome": ref_genome,
                "one_plasmid_file_path": gb , 
                "bowtie_path":"/home/yanghe/software/bowtie",   
                "no_ccdb_plasmid":"",  
                "no_sgRNA_plasmid":"",
                "plasmid_metod":plasmid_metod,
                "scene":scene,
                'sgRNA_result':{}, 

                "uha_dha_config": {
                    "max_right_arm_seq_length": 145,  
                    "max_left_arm_seq_length": 145,     
                    "min_left_arm_seq_length": 145,   
                    "min_right_arm_seq_length": 145     
                    },
                "plasmid_label":{
                            "ccdb_label":"ccdB",  
                            "promoter_terminator_label":"gRNA",
                            "n_20_label":"N20",
                            "promoter_label":"promoter"
                    },
                        
                "primer_json":{},
                "region_label":"", 
                "sgRNA_primer_json":{},
                "ccdb_primer_json":{},   
                "sgRNA_region_label":"",
                "ccdb_region_label":"",  

                        
                "enzyme":{
                            "enzyme_name":"BbsI",
                            "gap_sequence":"AA",    
                            "protection_sequence":"CCA"   
                        },      
                        
                "UHA_ARGS":{
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,  
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        }, 
                "SEQ_ALTERED_ARGS":{
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,  
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },
                "DHA_ARGS":{
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,
                            "PRIMER_MAX_TM": 75,
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },

                "PLASMID_Q_ARGS":{
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,  
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },    
                "GENOME_Q_ARGS":{
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,     
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },  
                "UP_SGRNA_ARGS": {
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,  
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },
                "DOWN_SGRNA_ARGS": {
                            "PRIMER_OPT_TM": 65,
                            "PRIMER_MIN_TM": 55,  
                            "PRIMER_MAX_TM": 75,    
                            "PRIMER_MIN_GC": 20,
                            'PRIMER_OPT_GC':65,
                            "PRIMER_MAX_GC": 80,
                            'PRIMER_MIN_SIZE':15,
                            'PRIMER_MAX_SIZE':25,
                            'PRIMER_OPT_SIZE':18, 
                        },
            }

    a = edit_main.main(data)
    print(a)

def excute_whole_process(process_id , opt_type, ref_genome_path, editInfo_input_csv, gb, editInfo_output_dir, chopchop_output_dir, edit_sequence_desgin_output_dir, analysis ):

    ref_genome = ref_genome_path
    editInfo_input_csv = editInfo_input_csv
    editInfo_output_dir = editInfo_output_dir
    scene = "both_sgRNA_primer"
    chopchop_input_csv = os.path.join(editInfo_output_dir,"info_input.csv" )
    chopchop_output_dir = chopchop_output_dir
    edit_sequence_desgin_input_csv = os.path.join(editInfo_output_dir,"info_input.csv" )
    sgRNA_csv = os.path.join(chopchop_output_dir,"sgRNA.csv" )
    edit_sequence_desgin_output_dir = edit_sequence_desgin_output_dir
    gb = gb
    plasmid_metod="0" 
  
    import time
    start_time_1_0 = time.time()  
    excute_data_preproccessing(editInfo_input_csv,ref_genome, editInfo_output_dir, scene)  
    end_time_1_0 = time.time() 
    exec_time_1_0 = end_time_1_0 - start_time_1_0
    print("函数的执行时间为：", exec_time_1_0, "秒")

    start_time_1_1 = time.time()  
    excute_chopchop(chopchop_input_csv, ref_genome, chopchop_output_dir, chopchop_config="")  
    end_time_1_1 = time.time() 
    exec_time_1_1 = end_time_1_1 - start_time_1_1
    print("函数的执行时间为：", exec_time_1_1, "秒")

    start_time_1_2 = time.time() 
    excute_edit_sequence_design( edit_sequence_desgin_input_csv, sgRNA_csv, edit_sequence_desgin_output_dir, ref_genome, gb, scene, plasmid_metod, analysis)
    end_time_1_2 = time.time() 
    exec_time_1_2 = end_time_1_2 - start_time_1_2
    print("函数的执行时间为：", exec_time_1_2, "秒")


    exec_time = end_time_1_2 - start_time_1_0
    time_dict = {
                    'strain':  basename(splitext(ref_genome_path)[0]),
                    'design target': process_id,
                    'preproccessing_time':exec_time_1_0,
                    'chopchop_time': exec_time_1_1,
                    'edit_sequence_design_time':exec_time_1_2,
                    'time': exec_time
                }

    exec_time = f'{exec_time_1_0},{exec_time_1_1},{exec_time_1_2},{exec_time}'
    
    pd.DataFrame([time_dict]).to_csv(  os.path.join(  dirname(ref_genome_path), ( basename(splitext(ref_genome_path)[0])  ),str(process_id), str(process_id)+f"_{opt_type}_time.csv"), index=False   )
    

    return exec_time

def execute_process(process_id, opt_type, ref_genome_path, editInfo_input_csv, gb, editInfo_output_dir, chopchop_output_dir, edit_sequence_desgin_output_dir,result_queue, analysis):

    print(f"Process {process_id} started")
    result = excute_whole_process(process_id , opt_type, ref_genome_path, editInfo_input_csv, gb, editInfo_output_dir, chopchop_output_dir, edit_sequence_desgin_output_dir,analysis)
    result_queue.put((process_id, result))
    print(f"Process {process_id} finished")

def oneprocessing_main(input_dir ,genome_name, whole_process_name, opt_type, index, analysis):

    #创建突变
    # input_dir = "/home/yanghe/program/data_preprocessing/genome_test"
    
    # genome_name_list = ['Escherichia_coli.fna', "Corynebacterium_glutamicum.fna","Pseudomonas_aeruginosa.fna","Bacillus_subtilis.fna","Salmonella_typhimurium.fna"]

    genome_name_list = ["Aspergillus_niger.fna"]   
    genome_dict_editInfo_csv = generate_batch_genome_mutation(genome_name_list, input_dir)
    pp(genome_dict_editInfo_csv)
    print("5个微生物基因组，50-500个编辑位点，100bp内的任意编辑类型，随机突变创建成功！")


    whole_process_dir = os.path.join(input_dir, genome_name, index, whole_process_name )

    ref_genome_path = os.path.join( input_dir,   genome_name+".fna" )
    editInfo_input_csv = os.path.join( input_dir, genome_name, index , index+".csv")
        
    editInfo_output_dir = os.path.join(whole_process_dir,"data_preprocessing")
    chopchop_output_dir = os.path.join(whole_process_dir,"chopchop")
    edit_sequence_desgin_output_dir = os.path.join(whole_process_dir, "edit_sequence_desgin")
  
    gb = "/home/yanghe/program/edit_sequence_design/input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori .gb"

    exec_time = excute_whole_process(index, opt_type, ref_genome_path, editInfo_input_csv, gb, editInfo_output_dir, chopchop_output_dir, edit_sequence_desgin_output_dir,analysis)

    preproccessing_time, chopchop_time, edit_sequence_design_time, all_time = exec_time.split(",")

    time_dict = {  
                    'strain':genome_name,
                    'design target': index,
                    'preproccessing_time':preproccessing_time,
                    'chopchop_time':chopchop_time,  
                    'edit_sequence_design_time':edit_sequence_design_time,
                    'all_time': all_time
                }
    
    time_df = pd.DataFrame( [time_dict] ) 

    time_df.to_csv(f"/home/yanghe/program/data_preprocessing/genome_test/{genome_name}.csv")



    return time_dict
   
def multiprocessing_main( input_dir, genome_name, whole_process_name, opt_type, genome_dict_editInfo_csv, analysis): 
    
    # 启动 50 个进程,执行全流程 
    result_queue = multiprocessing.Queue()
    processes = []   
    # genome_name = "Pseudomonas_aeruginosa"  


    ref_genome_path = os.path.join( input_dir,   genome_name+".fna" )
    gb = "/home/yanghe/program/edit_sequence_design/input/pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori .gb"


    for i,value in enumerate( genome_dict_editInfo_csv[genome_name] ):
    
        index = basename(value)

        editInfo_input_csv = os.path.join( input_dir, genome_name, index ,index+".csv")  
        whole_process_dir = os.path.join(input_dir, genome_name, index, whole_process_name)  
        editInfo_output_dir = os.path.join(whole_process_dir,"data_preprocessing")
        chopchop_output_dir = os.path.join(whole_process_dir,"chopchop")
        edit_sequence_desgin_output_dir = os.path.join(whole_process_dir, "edit_sequence_desgin")
       
        process = multiprocessing.Process(target=execute_process, args=(index, opt_type, ref_genome_path, editInfo_input_csv, gb, editInfo_output_dir, chopchop_output_dir, edit_sequence_desgin_output_dir, result_queue, analysis))
         
        processes.append(process)   
        process.start()
   
    # 等待所有进程完成
    for process in processes:
        process.join() 

    # 从队列中获取每个进程的结果
    results = []
    while not result_queue.empty():
        result_item = result_queue.get()
        design_target= result_item[0]
        preproccessing_time, chopchop_time, edit_sequence_design_time, all_time =   result_item[1].split(",")

        time_dict = {  
                    'strain':genome_name,
                    'design target': design_target,
                    'preproccessing_time':preproccessing_time,
                    'chopchop_time':chopchop_time,  
                    'edit_sequence_design_time':edit_sequence_design_time,
                    'all_time': all_time
                }
        results.append(time_dict)
    
    time_df = pd.DataFrame( results ) 

    # 打印每个进程的结果   
    pp( results )
    print("All processes have finished")
    
    return time_df
 
def patching_multiprocessing(input_dir,whole_process_name, genome_name, opt_type):

    # genome_name_list = ['Escherichia_coli.fna', "Corynebacterium_glutamicum.fna","Pseudomonas_aeruginosa.fna","Bacillus_subtilis.fna","Salmonella_typhimurium.fna"]
    genome_name_list = ['Yarrowia_lipolytica_CLIB_122.fna','Neurospora_crassa_OR74A.fna','Schizosaccharomyces_pombe_972h-.fna', 'Aspergillus_niger.fna']
    genome_dict_editInfo_csv = generate_batch_genome_mutation(genome_name_list, input_dir)
    pp(genome_dict_editInfo_csv)
    print("5个微生物基因组，40-400个编辑位点，100bp内的任意编辑类型，随机突变创建成功！")

    genome_dict_editInfo_csv[genome_name]

    new_genome_dict_editInfo_csv = {
                        genome_name:[]
    }

    for i,value in enumerate( genome_dict_editInfo_csv[genome_name] ):

        exists_path = os.path.join(value, f"{  basename(value) }_{opt_type}_time.csv" )

        
        if os.path.exists(exists_path):
            df = pd.read_csv(exists_path)
            if len(df.columns)>3:
                pass
            else:
                new_genome_dict_editInfo_csv[genome_name].append(value) 
        else:

            new_genome_dict_editInfo_csv[genome_name].append(value)  
    
    pp(new_genome_dict_editInfo_csv)
        
    time_df =  multiprocessing_main(input_dir, genome_name, whole_process_name, opt_type, new_genome_dict_editInfo_csv)   

    time_df.to_csv( os.path.join(input_dir,f"{genome_name}_{opt_type}.csv"), mode='a', header=False, index=False )

def execute_multiprocessing_main(input_dir, whole_process_name, genome_name, opt_type,analysis ):
    # genome_name_list = ['Escherichia_coli.fna', "Corynebacterium_glutamicum.fna","Pseudomonas_aeruginosa.fna","Bacillus_subtilis.fna","Salmonella_typhimurium.fna"]
    genome_name_list = ['Yarrowia_lipolytica_CLIB_122.fna'] 
    # genome_name_list = ['Aspergillus_niger.fna'] 
    genome_dict_editInfo_csv = generate_batch_genome_mutation(genome_name_list, input_dir)    
    pp(genome_dict_editInfo_csv)
    print("5个微生物基因组，40-400个编辑位点，100bp内的任意编辑类型，随机突变创建成功！")
    
    time_df = multiprocessing_main(input_dir, genome_name, whole_process_name, opt_type, genome_dict_editInfo_csv, analysis)
    time_df.to_csv( os.path.join(input_dir, f"{genome_name}_{opt_type}.csv"),index=False )   

def count_time(input_dir, genome_name,opt_type):

    # genome_name_list = ['Escherichia_coli.fna', "Corynebacterium_glutamicum.fna","Pseudomonas_aeruginosa.fna","Bacillus_subtilis.fna","Salmonella_typhimurium.fna"]
    genome_name_list = ['Yarrowia_lipolytica_CLIB_122.fna', 'Schizosaccharomyces_pombe_972h-.fna']
    # genome_name_list = ['Yarrowia_lipolytica_CLIB_122','Neurospora_crassa_OR74A','Schizosaccharomyces_pombe_972h-', 'Aspergillus_niger'] 
    # genome_name_list = ["Saccharomyces_cerevisiae_S288C.fna"]
    # genome_name_list = [ 'Aspergillus_niger.fna'] 
    genome_dict_editInfo_csv = generate_batch_genome_mutation(genome_name_list, input_dir)
    pp(genome_dict_editInfo_csv)
    print("5个微生物基因组，40-400个编辑位点，100bp内的任意编辑类型，随机突变创建成功！")  

    time_df = pd.DataFrame()

    for i,value in enumerate( genome_dict_editInfo_csv[genome_name] ):

        index = basename(value)
        
        time_path = os.path.join( value, f"{index}_before_time.csv")
        time_df = time_df.append(pd.read_csv(time_path)) 

    time_df.to_csv( os.path.join(input_dir,f"{genome_name}_{opt_type}.csv"),index=False)

def main():
      
    #创建突变  
    input_dir = "/home/yanghe/program/data_preprocessing/genome_test"  

    # genome_names = ['Aspergillus_niger']       
    # genome_names = ["Saccharomyces_cerevisiae_S288C"]   'Aspergillus_niger','Neurospora_crassa_OR74A'
    # genome_names = ["Escherichia_coli","Corynebacterium_glutamicum","Pseudomonas_aeruginosa","Bacillus_subtilis","Salmonella_typhimurium"]  'Saccharomyces_cerevisiae_S288C'

    # genome_name_list = ['Yarrowia_lipolytica_CLIB_122']   
    
    # for genome_name in genome_name_list:  
    #     opt_type = "before"  
    #     whole_process_name = "opt_before_whole_process_output"  
    #     analysis= False  # True：after，False：before
    #     execute_multiprocessing_main(input_dir, whole_process_name, genome_name, opt_type, analysis) 
    #     print(f"{genome_name},{opt_type}执行完毕")

    genome_name = "Aspergillus_niger"
    opt_type="before"  
    index = "450-1"
    whole_process_name = "opt_before_whole_process_output"
    analysis=True 
    oneprocessing_main( input_dir, genome_name, whole_process_name,  opt_type, index, analysis)


    genome_name = "Aspergillus_niger"
    opt_type="after"  
    index = "450-1"
    whole_process_name = "opt_after_whole_process_output"  
    analysis=False 
    oneprocessing_main( input_dir, genome_name, whole_process_name,  opt_type, index, analysis) 

        
        # opt_type = "after"
        # whole_process_name = "opt_after_whole_process_output"    
        # analysis= True  # True：after，False：before
        # execute_multiprocessing_main(input_dir, whole_process_name, genome_name, opt_type, analysis)  
        # print(f"{genome_name},{opt_type}执行完毕")  

        # patching_multiprocessing(input_dir,whole_process_name, genome_name, opt_type)      

        # count_time(input_dir, genome_name, opt_type)       

        #分析流程
        # import test
        # analysis_path =  f"/home/yanghe/program/test/analysis_{i}/"
        # test_analysis.main(analysis_path)

if __name__ == '__main__':
           
    main()





    