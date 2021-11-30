configfile: "config.yaml"

import io
import os
from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import numpy as np
import pathlib
import random
import pickle


rule salmon_index:
    input:
        raw_reads_1 = os.path.join(config["sample_dir"], "raw_files",
                                   "{sample_num}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["sample_dir"], "raw_files",
                                   "{sample_num}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["sample_dir"], "assembly_files",
                                     "{sample_num}.fasta")
    output:
        os.path.join(config["outputdir"], "01-salmon_index",\
                     "{sample_num}_index", "versionInfo.json")   
    params:
        libtype = "A",
        indexname = os.path.join(config["outputdir"], "01-salmon_index",\
                     "{sample_num}_index"),
        kval = 31
    conda:
        os.path.join("..", "envs", "salmon-env.yaml")
    log:
        err = os.path.join(config["output_dir"], "logs", "salmon", "{sample_num}_index.err"),
        out = os.path.join(config["output_dir"], "logs", "salmon", "{sample_num}_index.out")
    shell:
        '''
        salmon index -t {input.mock_assembly} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        '''

rule salmon_quant:
    input:
        raw_reads_1 = os.path.join(config["sample_dir"], "raw_files",
                                   "{sample_num}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["sample_dir"], "raw_files",
                                   "{sample_num}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["sample_dir"], "assembly_files",
                                     "{sample_num}.fasta"),
	index_complete = os.path.join(config["output_dir"], "01-salmon_index",\
                                      "{sample_num}_index", "versionInfo.json")  
    output: 
        os.path.join(config["output_dir"], "02-salmon_mapping",\
                     "{sample_num}_quant", "quant.sf")  
    params:
        libtype = "A",
        indexname = os.path.join(config["output_dir"], "02-salmon_mapping",\
                     "{sample_num}_index"),
        outdir = os.path.join(config["output_dir"], "02-salmon_mapping",\
                     "{sample_num}_quant")
    conda:
        os.path.join("..", "envs", "post-process-env.yaml")
    log:
        err = os.path.join(config["output_dir"], "logs", "salmon", "{sample_num}_quant.err"),
        out = os.path.join(config["output_dir"], "logs", "salmon", "{sample_num}_quant.out")
    shell:
        '''
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.raw_reads_1} -2 {input.raw_reads_2} --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        '''
