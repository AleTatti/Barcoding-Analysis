#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import subprocess
import datetime
from time import sleep
from Bio import SeqIO, SeqFeature, Entrez, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import MafftCommandline
from xml.dom import minidom
import re
import sys
from scipy.integrate import odeint
from statistics import mean, median, mode, stdev
import scipy.stats
from scipy import stats
import platform
import os
from os import listdir
from os.path import isfile, join
timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
import argparse


# In[2]:


#I define a function with an input the path of an excel in which I declare:
#the email for the NCBI use
#the path of the directory in which you want to work
#the executive file of blastn
#the executive file of mafft
#the executive file of iqtree
def read_configuration_file(filename):
    config = pd.read_excel(io=filename,
                      header = None,
                      names = ["TipologiaVariabile","ValoreVariabile"],
                          engine="openpyxl")
    email = config[config["TipologiaVariabile"]=="Email"]["ValoreVariabile"][0]
    path = config[config["TipologiaVariabile"]=="Path"]["ValoreVariabile"][1]
    blastn_exe = config[config["TipologiaVariabile"]=="blastn_exe"]["ValoreVariabile"][2]
    mafft_exe = config[config["TipologiaVariabile"]=="mafft_exe"]["ValoreVariabile"][3]
    iqtree = config[config["TipologiaVariabile"]=="iqtree"]["ValoreVariabile"][4]
    return [email, path, blastn_exe, mafft_exe, iqtree]      


# In[3]:


[Entrez.email, path, blastn_exe, mafft_exe, iqtree] = read_configuration_file("//Users//utente//projects//ConfigurationFile.xlsx")


# In[4]:


# I make a "TEMP" directory with a "Concatenate" directory inside it
dir = os.path.join(path,"TEMP","Concatenate")
if not os.path.exists(dir):
    os.mkdir(dir)


# In[5]:


#I define a function to take the first word of a string, based on ":"
def first_word(x):
    output= x.split(":")[0]
    return output
#I define a function to convert a string to list
def tolist(x):
    y = list(x.split(";"))
    y = [float(i) for i in y]
    return y


# In[6]:


# I define a function to remove a list of common geographical words
def EliminaParole(x):
    f = open(path + "Words.txt", "r") #https://github.com/AleTatti/Barcoding-Analysis/blob/main/words.txt
    Words = f.read()
    f.close
    stopwords = Words.lower().split("\n")
    parole = x.split(' ')
    paroleNew = [a for a in parole if a not in stopwords]
    paroleNew = [a for a in paroleNew if a not in ""]  
    return paroleNew


# In[7]:


class Create_Dataset(object):
    def __init__(self, path, blastn_exe, mafft_exe, iqtree): 
        self.path = path
        self.blastn_exe = blastn_exe
        self.mafft_exe = mafft_exe
        self.iqtree = iqtree
    
    def read_file(self,filename,tipologia): 
        if filename.startswith("//") or filename.startswith("\\"):
            f = open(filename, "r")
        else:
            f = open(self.path + filename, "r")
            if tipologia == "list":                                   
                self.filename_taxa = filename
            elif tipologia == "fasta" or tipologia == "aln":         
                self.filename_fasta = filename
        
        self.data_in_string = f.read()
        f.close()
        
    def transform_file_in_filter(self, filename, gene, rank_type):
        gene = gene.split(",")
        lista_gene = []
        for j in range(len(gene)):
            stringa_gene = ""
            gene_split = gene[j].split(" ")
            if gene[j] == "large subunit ribosomal RNA" or gene[j] == "LSU" or gene[j] == "28S" or gene[j]=="18S" or gene[j] =="SSU" or gene[j]=="small subunit ribosomal RNA":
                for i in range(len(gene_split)):
                    stringa_gene = stringa_gene + " AND " + gene_split[i] + "[Title]"  
                stringa_gene = stringa_gene + " NOT internal[Title]" + " NOT ITS2[Title]"
            else:
                for i in range(len(gene_split)):
                    stringa_gene = stringa_gene + " AND " + gene_split[i] + "[Title]"  
            lista_gene.append(stringa_gene)
            
        if filename.endswith(".txt"): 
            self.read_file(filename, tipologia = "list")
            dati_splittati = self.data_in_string.split("\n")
        else:
            dati_splittati = filename.split(",")
        self.dati_splittati = dati_splittati
        
        lista_taxa = []
        if rank_type == "A":
            for i in range(len(dati_splittati)):    
                lista_taxa.append(dati_splittati[i].split(" ")[0]) 
        elif rank_type == "S":
            for i in range(len(dati_splittati)):    
                lista_taxa.append(dati_splittati[i]) 
        else:
            print("specify the rank type: S = species; A = above species level")
        
        lista_taxa = list(set(lista_taxa))
        lista_filtri = []
        
        for j in range(len(lista_gene)):
            for i in range(len(lista_taxa)):
                lista_filtri.append("(\"" + lista_taxa[i] + "\"[Organism])" + lista_gene[j])
    
        self.lista_filtri = lista_filtri
        self.gene = gene  
        
        
    def find_accessid(self, par_retmax):
    
        lista_filtri = self.lista_filtri
        access_id = []
    
        for filtro in range(len(lista_filtri)):
            handle = Entrez.esearch(
              db="nucleotide", 
              term = lista_filtri[filtro], 
              idtype="acc", 
              retmax = par_retmax, 
              field="title")  

            record = Entrez.read(handle) 
            accession = record["IdList"]

            if(len(accession) == par_retmax):
                print('par_retmax is too low, there are more sequences to find') 

            print(lista_filtri[filtro])
            print((len(accession)), " sequences found")

            access_id = access_id + accession
              
        self.accessid_from_keyword = list(set(access_id))
        print(len(self.accessid_from_keyword), "Accession Numbers found from the keyword search")
        
    def keyword_search(self, filename, gene, rank_type):
        print('searching by keyword')
        self.transform_file_in_filter(filename,gene,rank_type)
        self.find_accessid(par_retmax=10000)
         
    def format_fasta(self,filename): 
        if filename.endswith("aln"):
            
            self.read_file(filename,tipologia="aln") 
        else:
            self.read_file(filename,tipologia="fasta") 
        DatiSplittati = self.data_in_string.split("\n") 
        stringa = ""
        for i in range(len(DatiSplittati)):
                if(len(DatiSplittati[i])>0): 
                    if(DatiSplittati[i][0] == ">"): 
                        stringa = stringa + "\n" + DatiSplittati[i] + "\n"
                if(len(DatiSplittati[i])>0): 
                    if(DatiSplittati[i][0] != ">"): 
                        stringa = stringa + DatiSplittati[i] 
                if(len(DatiSplittati[i])==0): 
                        stringa = stringa
        self.data_in_string = stringa[1:]
        
    def print_file_fasta(self):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        dir = os.path.join(path,"TEMP")
        if not os.path.exists(dir):
            os.mkdir(dir)
     
        if platform.system() == "Darwin" or platform.system() == "Linux":
            if self.filename_fasta.startswith("//") or self.filename_fasta.startswith("TEMP//"): 
                self.filename_fasta = self.filename_fasta.split("//")[-1] 
            file_output = self.path + "TEMP//" + self.filename_fasta + "_" + timestamp + ".fasta"
                
        if platform.system() == "Windows":
            if self.filename_fasta.startswith("\\") or self.filename_fasta.startswith("TEMP\\"): 
                self.filename_fasta = self.filename_fasta.split("\\")[-1] 
            file_output = self.path + "TEMP\\" + self.filename_fasta + "_" + timestamp + ".fasta"
        
        self.filefasta_output = file_output

        output = open(file_output, "w")
        output.write(self.data_in_string)
        output.close()

    
    def run_blast(self,query):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + '_blasting the input sequence/sequences...')
        
        if str(query).endswith(".fasta"):
            self.format_fasta(query)
            self.print_file_fasta()

            query = self.filefasta_output

        else:
            query = query
         
        cmd =  self.blastn_exe
        out =  str(query) + timestamp + ".csv"

        self.output_post_blast = out
        #outfmt: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        cline = NcbiblastnCommandline(remote = True, query = query, db="nt", strand="plus",
                                      evalue=0.001, out=out, outfmt=6, cmd = cmd)

        cline()
        
    def read_blast(self,threshold,csv="none"):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + "_reading the csv file obtained from blast...")
        
        if (csv == "none"):
            FromBlastDf = pd.read_csv(self.output_post_blast,
                          sep = "\t",
                          header=None)
        else:
            FromBlastDf = pd.read_csv(self.path+csv,
                          sep = ",",
                          header=None)

        FromBlastDf = FromBlastDf[FromBlastDf[2]>threshold] 
        self.accessid_from_blast = list(set(list(FromBlastDf[1])))
        self.FromBlastDf = FromBlastDf
        
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + "_now you can read accessid_from_blast")
        print(len(self.accessid_from_blast), "Accession Numbers found from the blast search")
          
    def execute_blast(self,query,threshold,csv):
        self.run_blast(query)
        print(self.path)
        self.read_blast(threshold,csv)
    
    
    def merge_accessid(self):
        self.accessid_final = list(set(self.accessid_from_blast + self.accessid_from_keyword))
        print("A total of", len(self.accessid_final), "Accession Numbers found")
        
        
    def find_from_accessid(self, query, threshold, csv, filename, id_type, gene, rank_type):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        if id_type == 'total': 
            self.execute_blast(query,threshold,csv)
            self.keyword_search(filename,gene,rank_type)
            self.merge_accessid()
            if len(self.accessid_final)==0:
                print("no sequences found")
                pass
            else:
                print(timestamp, " searching for the access ids from blast and keyword search on the GenBank database")
                handle = Entrez.efetch(db="nucleotide", 
                                           id= self.accessid_final, 
                                           rettype="gb", 
                                           retmode="xml") 
        elif id_type == 'K': 
            self.keyword_search(filename,gene,rank_type)
            if len(self.accessid_from_keyword)==0:
                print("no sequences found")
                pass
            else:
                print(timestamp, " searching for the access ids from the keyword search on the GenBank database")
                handle = Entrez.efetch(db="nucleotide", 
                                           id= self.accessid_from_keyword,
                                           rettype="gb", 
                                           retmode="xml") 
            
        elif id_type == 'B': 
            self.execute_blast(query,threshold,csv)
            print(timestamp, " searching for the access ids from the blast search on the GenBank database")
            handle = Entrez.efetch(db="nucleotide", 
                                       id= self.accessid_from_blast,
                                       rettype="gb", 
                                       retmode="xml") 
                        
        elif id_type == 'aln': 
            print(timestamp, " reading the alignment file")
            print(timestamp, " searching for the access ids from the fasta file on the GenBank database")
            handle = Entrez.efetch(db="nucleotide", 
                                   id=self.df_gene['accessid'], 
                                   rettype="gb", 
                                   retmode="xml") 
            
        elif id_type == 'file': 
            print(timestamp + " reading the access id file")
            self.read_file(filename,tipologia ="list")
            access_id_list = self.data_in_string.split("\n")
            print(timestamp + " searching for the access ids of the input file on the GenBank database")
            handle = Entrez.efetch(db="nucleotide", 
                                       id = access_id_list,
                                       rettype="gb", 
                                       retmode="xml")  
            
        elif id_type == 'list':
            print("searching for the access id list on the GenBank database")
            handle = Entrez.efetch(db="nucleotide", 
                                   id = filename, 
                                   rettype="gb", 
                                   retmode="xml") 
        if id_type == "K":
            if len(self.accessid_from_keyword)==0:
                print("no sequences found")
                pass
            else:
                record = Entrez.read(handle)
                handle.close()
                self.record = record
        if id_type == "total":
            if len(self.accessid_final)==0:
                print("no sequences found")
                pass
            else:
                record = Entrez.read(handle)
                handle.close()
                self.record = record

        lista = []
        for E in range(len(record)):
            if(('GBSeq_sequence' in record[E]) == True): 
                accessid = record[E]["GBSeq_accession-version"]
                ElencoFeatures = record[E]["GBSeq_feature-table"][0]["GBFeature_quals"]
                species = ""
                type_material = ""
                country = ""
                collection_date = ""
                strain = ""
                cultivar = ""
                isolate = ""
                clone = ""
                culture_collection = ""
                specimen_voucher = ""
                isolation_source = ""
                host = ""
                db_xref = ""
                lat_lon = ""
                altitude = ""
                collected_by = ""
                identified_by = ""
                note = ""
                PCR_primers = ""
                searching_gene = record[E]["GBSeq_feature-table"]
                self.searching_gene = searching_gene
                sequence = record[E]["GBSeq_sequence"]
                for y in range(len(searching_gene)): 
                    if(('GBFeature_quals' in searching_gene[y]) == True):
                        for z in range(len(searching_gene[y]["GBFeature_quals"])): 
                            try:
                                if(searching_gene[y]["GBFeature_quals"][z]['GBQualifier_value'] in gene):
                                    if searching_gene[y]["GBFeature_quals"][z]['GBQualifier_name'] == "transcription":
                                        sequence = searching_gene[y]["GBFeature_quals"][z]["GBQualifier_value"]
                                        break
                                    start = int(record[E]["GBSeq_feature-table"][y]['GBFeature_location'].replace("<","").replace(">","").replace("complement(","").replace("(","").split("..")[0])
                                    end = int(record[E]["GBSeq_feature-table"][y]['GBFeature_location'].replace("<","").replace(">","").replace("complement(","").replace(")","").split("..")[1])
                                    sequence = record[E]["GBSeq_sequence"][start:end]
                                    break
                            except:
                                pass
                         
                            
                for i in range(len(ElencoFeatures)):
                        if( ElencoFeatures[i]["GBQualifier_name"] == "organism"):
                            species = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "country"):
                            country = ElencoFeatures[i]["GBQualifier_value"]    
                        if( ElencoFeatures[i]["GBQualifier_name"] == "lat_lon"):
                            lat_lon = ElencoFeatures[i]["GBQualifier_value"] 
                        if(ElencoFeatures[i]["GBQualifier_name"] == 'altitude'):
                            altitude = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "collection_date"):   
                            collection_date =  ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "isolation_source"):
                            isolation_source = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "host"):
                            host = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "cultivar"):
                            cultivar = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "isolate"):
                            isolate = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "strain"):
                            strain = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "clone"):
                            clone = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "culture_collection"):
                            culture_collection = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "specimen_voucher"):
                            specimen_voucher = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "db_xref"):   
                            db_xref =  ElencoFeatures[i]["GBQualifier_value"] 
                        if( ElencoFeatures[i]["GBQualifier_name"] == "collected_by"):   
                            collected_by =  ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "identified_by"):
                            identified_by = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "type_material"):
                            type_material = ElencoFeatures[i]["GBQualifier_value"]
                        if( ElencoFeatures[i]["GBQualifier_name"] == "note"):
                            note = ElencoFeatures[i]["GBQualifier_value"]
                        if ( ElencoFeatures[i]["GBQualifier_name"] == "PCR_primers"):
                            PCR_primers = ElencoFeatures[i]["GBQualifier_value"]
                riga = [accessid, species, type_material, country, lat_lon, altitude, isolation_source, host, 
                        db_xref, isolate, cultivar, strain, clone, culture_collection, specimen_voucher, 
                        collected_by, identified_by, collection_date, PCR_primers, note, sequence]
                lista.append(riga)
            else:
                continue
                  
        df = pd.DataFrame(lista, columns=["accessid", "species","type_material","country", "lat_lon", "altitude","isolation_source", "host", 
                                          "db_xref", "isolate", "cultivar","strain", "clone", "culture_collection","specimen_voucher", 
                                          "collected_by", "identified_by","collection_date", "PCR_primers", "note","sequence"])
        df["cleaned_species"] = ""
        notsp = ["aff.","cf.","nr.","uncultured","UNVERIFIED", "Uncultured"]
        for i in range(len(df["cleaned_species"])):
             df["cleaned_species"][i] = df["species"][i].split("-")[0]
        df["cleaned_species"] = df["cleaned_species"].str.replace('[^a-zA-Z. ]','')
        for i in range(len(df["cleaned_species"])):
            for word in notsp:
                if (word in df["cleaned_species"][i].split(" ")):
                    df["cleaned_species"][i] = ' '.join(map(str,[a for a in df["species"][i].split(" ") if a != word]))
            if "sp." in df["cleaned_species"][i].split(" "):
                if len(df["cleaned_species"][i].split(" "))>2:
                    if len(df["cleaned_species"][i].split(" ")[2])<3 or df["cleaned_species"][i].split(" ")[2][0].isupper():
                        df["cleaned_species"][i] = df["cleaned_species"][i].split(" ")[0] + " sp."
                    else:
                        df["cleaned_species"][i] = df["cleaned_species"][i].split(" ")[0] + " " +  df["cleaned_species"][i].split(" ")[2]
            if (len(df["cleaned_species"][i].split(" ")) == 1) or (df["cleaned_species"][i].split(" ")[1][0].isupper()):
                df["cleaned_species"][i] = df["cleaned_species"][i].split(" ")[0] + " sp."
            if (len(df["cleaned_species"][i].split(" ")) > 2):
                df["cleaned_species"][i] = df["cleaned_species"][i].split(" ")[0] + " " + df["cleaned_species"][i].split(" ")[1]
        df["state"] = df['country'].apply(first_word)
        df = df.sort_values(by=['cleaned_species'])
        self.df = df
        print('''Dataframe "df" created with ''', len(df), " sequences")
                
    def Save_fasta(self, cleaning_level, n, G, Dataframe, out, header="default"):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, " saving the dataset fasta file..")
        Dataframe = Dataframe
        if str(Dataframe).endswith(".xlsx"):
            self.Open_Dataframe(filename=Dataframe, concat="No")
            Dataframe = self.MyDf
        if header == "default":
            if "concatseq" in Dataframe.columns:
                header = "accessid,Species"
            else:
                header = "accessid,species,state"
        header=header.split(",")
        self.header = header
        stringa = ""
        Dataframe = Dataframe.reset_index(drop = True)
        Dataframe['len'] = '' 
        if "aligned_sequence" in Dataframe.columns:
            Dataframe['sequence'] = Dataframe["aligned_sequence"]
        Dataframe0 = Dataframe.copy()
        if cleaning_level > 0:   
            print("Cleaning the dataset")
            print("Dataset starts with ", len(Dataframe), " sequences")
            if "concatseq" in Dataframe.columns:
                Dataframe["len"]=Dataframe["concatseq"].str.len()
                print("removing gap-only columns")
                lista_aligned_sequence = Dataframe["concatseq"].tolist()
                lista_di_lista_aligned_sequence=[]
                for i in range(len(lista_aligned_sequence)):
                    lista_di_lista_aligned_sequence.append(list(lista_aligned_sequence[i]))
                arr = np.array(lista_di_lista_aligned_sequence)
                posizioni_da_tenere = []
                for i in range(np.size(arr[0,:])):
                    if (np.char.count(arr[:,i],"-").sum() != np.size(arr[:,i])):
                        posizioni_da_tenere.append(i)
                arr_pulito = arr[:,posizioni_da_tenere]
                lista_pulita = []
                for i in range(len(arr_pulito)):
                    lista_pulita.append(''.join(list(arr_pulito[i])))
                Dataframe["concatseq"] = lista_pulita        
            else:
                Dataframe["len"]=Dataframe["sequence"].str.len()
            Dataframe = Dataframe.reset_index(drop=True)
            if len(Dataframe) > 1:
                print("removing too long sequences..")
                Dataframe
                for i in range(len(Dataframe)):
                    if Dataframe['len'][i] > (mean(Dataframe['len']) + 10*stdev(Dataframe['len'])):
                        print(str(Dataframe.loc[i,"accessid"]), " is too long. It could be a different set of primers, a whole genome or chromosome and will be excluded")
                        self.Dataframe = Dataframe
                        Dataframe = Dataframe[Dataframe['len'] < (mean(Dataframe['len']) + 10*stdev(Dataframe['len']))]  
            Dataframe = Dataframe.reset_index(drop=True)
            Dataframe1 = Dataframe.copy()
            print(len(Dataframe0)-len(Dataframe), "sequences deleted")
            if cleaning_level >= 2:                 
                print('deleting sequences with more than', n, ' of ambiguous nucleotides')
                Dataframe_c = Dataframe.copy
                if "concatseq" in Dataframe.columns:
                    Dataframe['concatseqN'] = Dataframe['concatseq'].str.replace('[^atcgATCG-]','N')
                    Dataframe['n_count'] = Dataframe['concatseqN'].str.count('N')
                else:
                    Dataframe['sequenceN'] = Dataframe['sequence'].str.replace('[^atcgATCG-]','N')
                    Dataframe['n_count'] = Dataframe['sequenceN'].str.count('N')
                
                Dataframe['x_n'] = Dataframe['n_count']/Dataframe['len']
                self.Dataframe1 = Dataframe.copy()
                Dataframe = Dataframe[Dataframe['x_n'] < n] 
                Dataframe = Dataframe.reset_index(drop=True)
                print(len(Dataframe1)-len(Dataframe), "sequences deleted")
                Dataframe2 = Dataframe.copy()
                if cleaning_level >= 3:
                    if "concatseq" not in Dataframe.columns: 
                        if "aligned_sequence" not in Dataframe.columns:
                            if len(Dataframe) > 1:
                                print("deleting sequences shorter than the mean - stdev of 33% of other sequences")
                                Dataframe = Dataframe[Dataframe['len'] > mean((Dataframe['len']*0.33)) - stdev(Dataframe['len']*0.33)]
                                Dataframe = Dataframe.reset_index(drop=True)
                                print(len(Dataframe2)-(len(Dataframe)),"sequences deleted")
                            else:
                                print("you need at least two sequences")
                        else:                          
                            print("trimming gap tail and head of the alignment")
                            Dataframe["gap_head"] = ''
                            Dataframe["gap_tail"] = ''
                            for i in range(len(Dataframe)):
                                head = re.search('\A-+', Dataframe["aligned_sequence"][i]) 
                                if head == None:
                                    Dataframe.loc[i, "gap_head"] = 0
                                else:
                                    Dataframe.loc[i, "gap_head"] = len(head.group(0))
                                tail = re.search('-+.$', Dataframe["aligned_sequence"][i]) 
                                if tail == None:
                                    Dataframe.loc[i,"gap_tail"] = 0
                                else:
                                    Dataframe.loc[i,"gap_tail"] = len(tail.group(0))
                            head_cut = Dataframe["gap_head"].value_counts()
                            self.DataframeN = Dataframe.copy()
                            head_trim = head_cut.index[0] 
                            if len(head_cut) > 1:
                                differenza = abs(head_cut.index[0] - head_cut.index[1])
                                if differenza < 5:
                                    print("you have two very similar head_gaps", head_cut.index[0], head_cut.index[1], ", the max number of gaps is chosen")
                                    head_trim = max(head_cut.index[0], head_cut.index[1])
                            tail_cut = Dataframe["gap_tail"].value_counts()
                            tail_trim = tail_cut.index[0] 
                            if(len(tail_cut)) > 1:
                                differenza = abs(tail_cut.index[0] - tail_cut.index[1])
                                if differenza < 5:
                                    print("you have two very similar tail_gaps", tail_cut.index[0], tail_cut.index[1], ", the max number of gaps is chosen")
                                    tail_trim = max(tail_cut.index[0], tail_cut.index[1])
                            Dataframe["trimmed_sequence"] = ''
                            for i in range(len(Dataframe)):
                                if tail_trim == 0:
                                    Dataframe.loc[i,"trimmed_sequence"] = Dataframe.loc[i,"aligned_sequence"][head_trim:]
                                else:
                                    Dataframe.loc[i,"trimmed_sequence"] = Dataframe.loc[i,"aligned_sequence"][head_trim:-tail_trim:]
                            print(head_trim, " positions removed from 5'")
                            print(tail_trim, " positions removed from 3'")       
                            Dataframe_preOdd = Dataframe
                            print("removing the odd one out")
                            lista_aligned_sequence = Dataframe["trimmed_sequence"].tolist()
                            lista_di_lista_aligned_sequence=[]
                            for i in range(len(lista_aligned_sequence)):
                                lista_di_lista_aligned_sequence.append(list(lista_aligned_sequence[i]))
                            arr = np.array(lista_di_lista_aligned_sequence)
                            posizione_access_id_da_escludere = []
                            contatore = 0
                            candidato_da_escludere = -1
                            for i in range(np.size(arr[0,:])):
                                if (np.char.count(arr[:,i],"-").sum() == np.size(arr[:,i])-1):
                                    candidato = int(np.argwhere(np.char.count(arr[:,i],"-")==0))
                                    if(candidato == candidato_da_escludere):
                                        contatore = contatore + 1
                                        if(contatore == 5):
                                            posizione_access_id_da_escludere.append(candidato)
                                    else:
                                        candidato_da_escludere = candidato
                                        contatore = 1
                                else:
                                    contatore = 0
                                    candidato_da_escludere = -1 
                            da_escludere = list(Dataframe.iloc[list(set(posizione_access_id_da_escludere))]["accessid"])        
                            Dataframe = Dataframe[~Dataframe["accessid"].isin(da_escludere)].reset_index(drop=True)
                            print(len(Dataframe_preOdd)-len(Dataframe), "sequences removed")
                            print("removing sequences with more gap than 80% of the length")
                            Dataframe["Gap_trim"] = Dataframe["trimmed_sequence"].str.count("-")
                            Dataframe["len"]=Dataframe["trimmed_sequence"].str.len()
                            #Dataframe = Dataframe[Dataframe['Gap_trim'] < (Dataframe['len']*0.8)] 
                            Dataframe = Dataframe.reset_index(drop=True)
                            self.Dataframe5 = Dataframe.copy()
                            print(len(Dataframe2)-len(Dataframe), "sequences deleted")                            
                            Dataframe3 = Dataframe.copy()
                            self.Dataframe3 = Dataframe.copy()                           
                            if cleaning_level == 4:
                                print("deleting all the sequences with more than" ,G, " gaps")
                                Dataframe["Gap_trim"] = Dataframe["trimmed_sequence"].str.count("-")
                                Dataframe = Dataframe[Dataframe['Gap_trim'] < Dataframe['len']*G] 
                                Dataframe = Dataframe.reset_index(drop=True)
                                print(len(Dataframe3)-len(Dataframe), "sequences deleted")
                                self.Dataframe6 = Dataframe.copy()
                                
                                print("excluding odd sequences")
                                lista_aligned_sequence = Dataframe["trimmed_sequence"].tolist()
                                lista_di_lista_aligned_sequence=[]
                                for n in range(len(Dataframe)):
                                    for i in range(len(lista_aligned_sequence)):
                                        lista_di_lista_aligned_sequence.append(list(lista_aligned_sequence[i]))
                                    arr = np.array(lista_di_lista_aligned_sequence)
                                    posizione_access_id_da_escludere = []
                                    contatore = 0
                                    candidato_da_escludere = -1
                                    for i in range(np.size(arr[0,:])):
                                        if (np.char.count(arr[:,i],"-").sum() == np.size(arr[:,i])-1):
                                            candidato = int(np.argwhere(np.char.count(arr[:,i],"-")==0))
                                            if(candidato == candidato_da_escludere):
                                                contatore = contatore + 1
                                                if(contatore == 3):
                                                    posizione_access_id_da_escludere.append(candidato)
                                            else:
                                                candidato_da_escludere = candidato
                                                contatore = 1
                                        else:
                                            contatore = 0
                                            candidato_da_escludere = -1 
                                    da_escludere = list(Dataframe.iloc[list(set(posizione_access_id_da_escludere))]["accessid"])        
                                    lunghezza_df = len(Dataframe)
                                    Dataframe = Dataframe[~Dataframe["accessid"].isin(da_escludere)].reset_index(drop=True)
                                    if(len(Dataframe)==lunghezza_df):
                                        break
                                print(len(self.Dataframe6)-len(Dataframe), "sequences deleted")   
                            Dataframe = Dataframe.reset_index(drop=True)                   
                            print("removing gap-only columns")
                            lista_aligned_sequence = Dataframe["trimmed_sequence"].tolist()
                            lista_di_lista_aligned_sequence=[]
                            for i in range(len(lista_aligned_sequence)):
                                lista_di_lista_aligned_sequence.append(list(lista_aligned_sequence[i]))
                            arr = np.array(lista_di_lista_aligned_sequence)
                            posizioni_da_tenere = []
                            for i in range(np.size(arr[0,:])):
                                if (np.char.count(arr[:,i],"-").sum() != np.size(arr[:,i])):
                                    posizioni_da_tenere.append(i)
                            arr_pulito = arr[:,posizioni_da_tenere]
                            lista_pulita = []
                            for i in range(len(arr_pulito)):
                                lista_pulita.append(''.join(list(arr_pulito[i])))
                            Dataframe["trimmed_sequence"] = lista_pulita            
                    if "concatseq" in Dataframe.columns: 
                        print("deleting sequences with more than 60% of gaps")
                        Dataframe = Dataframe[Dataframe['Gap'] < (Dataframe['len']*0.6)] 
                        print(len(Dataframe2)-len(Dataframe), "sequences deleted")
                        Dataframe3 = Dataframe.copy()
                        Dataframe = Dataframe.reset_index(drop=True)  
                    print("a total of ", len(Dataframe0)-len(Dataframe), "sequences deleted")
                    
                    
                    
                    
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, "saving the formatted fasta file")    
        for i in range(len(Dataframe)):
            stringa = stringa + ">" 
            for j in range(len(header)):
                stringa = stringa + str(Dataframe.loc[i,header[j]]) + "_"     
            stringa = re.sub('_+','_',stringa)
            stringa = stringa.strip("_")
            if "concatseq" not in Dataframe.columns: 
                if "aligned_sequence" in Dataframe.columns:
                    if cleaning_level >= 3:
                        stringa = stringa + "\n" + Dataframe.loc[i,"trimmed_sequence"] + "\n"
                    else:
                        stringa = stringa + "\n" + Dataframe.loc[i,"aligned_sequence"] + "\n"
                else:
                    stringa = stringa + "\n" + Dataframe.loc[i,"sequence"] + "\n"
            elif "concatseq" in Dataframe.columns: 
                if cleaning_level >= 3:
                    stringa = stringa + "\n" + Dataframe.loc[i,"trimmed_sequence"] + "\n"
                else:
                    stringa = stringa + "\n" + Dataframe.loc[i,"concatseq"] + "\n"
            stringa = re.sub('[^>a-zA-Z0-9.r\n-]','_',stringa)       
        stringa = re.sub('_+','_',stringa)
        self.file_fasta = stringa
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        file_output = path + out + "_" + timestamp + "_Dataset.fasta"
        self.fasta_output = file_output
        o = open(file_output, "w")
        o.write(stringa)
        o.close
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, " the fasta file is ready")
        

    def Build_Dataset(self, query, threshold, csv, filename, id_type, gene, rank_type = "S", cleaning_level=1, n=0.01, G=0.33):
            if len(gene.split(";"))>1:
                if id_type == "total":
                    print("Blast search plus keyword search need only one gene")
                else:
                    for i in range(len(gene.split(";"))):
                        self.find_from_accessid(query, threshold, csv, filename, id_type, gene.split(";")[i], rank_type)
                        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                        print(timestamp, "saving the dataset in a excel dataframe")
                        self.df.to_excel(path + "Dataframe_" + gene.split(";")[i].split(",")[0].replace(" ","_")+"_"+timestamp + ".xlsx",index = False)
                        df1 = self.df.copy()
                        out = filename.split(" ")[0] + "_" + gene.split(";")[i].split(",")[0].replace(" ","_")
                        self.Save_fasta(cleaning_level,n,G,Dataframe=df1,out=out,header="default")
            else:
                self.find_from_accessid(query, threshold, csv, filename, id_type, gene, rank_type)
                timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                print(timestamp, "saving the dataset in a excel dataframe")
                self.df.to_excel(path + "Dataframe_" + gene.split(",")[0].replace(" ","_") + "_" + timestamp + ".xlsx",index = False)
                df1 = self.df.copy()
                out = "Dataset_" + filename.split(" ")[0] + "_" + gene.split(",")[0].replace(" ","_")
                self.Save_fasta(cleaning_level,n,G,Dataframe=df1,out=out,header="default")

    def Align_fasta(self,filename):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, '''"aligning the dataset with mafft with the "adjustdirectionaccurately" option.."''')
        in_file = self.path + filename
        mafft_cline = MafftCommandline(mafft_exe, input=in_file, adjustdirectionaccurately= "on")
        print(mafft_cline)
        stdout, stderr = mafft_cline()
        with open(path + filename, "w") as handle:
            handle.write(stdout)
        align = AlignIO.read(path + filename, "fasta")
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, "the alignment is ready")
            
    def Open_Dataframe(self, filename, concat):      
        if concat == "Y":
            if platform.system() == "Darwin" or platform.system() == "Linux":
                MyDf = pd.read_excel(self.path + "Concatenate//" + filename, index_col=None, header = 0, keep_default_na=False, engine="openpyxl")
            if platform.system() == "Windows":
                MyDf = pd.read_excel(self.path + "Concatenate\\" + filename, index_col=None, header = 0, keep_default_na=False, engine="openpyxl")
        else:
            MyDf = pd.read_excel(self.path + filename, index_col=None, 
                                 header = 0, keep_default_na=False, engine="openpyxl")
        self.MyDf = MyDf  
        #print("The data frame is available as MyDf")
        
    def read_aln(self,alignment_file):
        filename = alignment_file
        self.format_fasta(filename)
        self.df_gene = pd.DataFrame(columns = ['accessid','aligned_sequence'])
        splitted_data = self.data_in_string.split(">")
        splitted_data.pop(0)
        data  = []
        for i in range(len(splitted_data)):
            header = splitted_data[i].split("\n")[0]
            prefix = header.split("_")[0]
            seq = splitted_data[i].split("\n")[1]
            if(len(prefix) == 0): 
                if(len(header.split("_")[2]) < 3):
                        data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                        'aligned_sequence': seq}]
                else:
                        data = data + [{'accessid': header.split("_")[2], 
                                        'aligned_sequence': seq}]
            else:
                if(len(prefix) < 3):
                    data = data + [{'accessid': prefix + "_" + header.split("_")[1],
                                    'aligned_sequence': seq}]
                else:
                    data = data + [{'accessid': prefix,
                                    'aligned_sequence': seq}]                                
        self.df_gene = self.df_gene.append(data,ignore_index=True)
        
        
    def merge_df_to_aln(self, df_type):
        if df_type == 'df':
            if(len(self.df) < len(self.df_gene)):
                print("have you deleted some sequences in the alignment? in that case never mind; otherwise check it")
            df = pd.merge(self.df, self.df_gene, how = 'inner',left_on = "accessid", 
                          right_on = "accessid",suffixes=('', ''))
            self.df = df
        elif df_type == 'MyDf':
            if(len(self.MyDf) < len(self.df_gene)):
                print("have you deleted some sequences in the alignment? in that case never mind; otherwise check it")
            MyDf = pd.merge(self.MyDf, self.df_gene, how = 'inner',left_on = "accessid", 
                            right_on = "accessid",suffixes=('', ''))
            self.df = MyDf
        else:
            print("ERROR: specify df meaning a new created Data Frame or MyDf an uploaded one")
              
    def correcting_info(self, organism, df_number):
        if organism == "Wolbachia":
            for i in range(len(self.df)):
                if ((self.df['species'][i].split(" ")[1] == 'pipientis') 
                    or (self.df['species'][i].split(" ")[1] == 'massiliensis') 
                    or (self.df['species'][i].split(" ")[1] == 'sp.') 
                    or (self.df['species'][i].split(" ")[2] == 'sp.')):
                    if (self.df['host'][i]) != '':
                        if len(self.df['host'][i].split(" ")) > 1:
                            self.df['species'][i] = 'Wolbachia endosymbiont of ' 
                            + self.df['host'][i].split(" ")[0] + " " + self.df['host'][i].split(" ")[1]
                        else:
                            self.df['species'][i] = 'Wolbachia endosymbiont of ' 
                            + self.df['host'][i].split(" ")[0]
                    else:
                        self.df['species'][i] = 'Wolbachia sp.'
            self.df['Species'] = self.df['species'].str.replace(' ','_') 
                   
        if organism == 'Fungi':
            bad_words = ["INATURALIST", "OBSERVER", "MUSHROOMOBSERVER", "HERBARIUM", "OBSERVATIONS"]
            not_too_bad_words = ["PBM", "UBC", "STRAIN", "STRAING", "CBS", "ISOLATE", "CLONE"]
            columns = ['specimen_voucher','strain', 'isolate',"clone" ]
            for col in columns:
                for i in range(len(self.df[col])):
                    self.df[col][i] = re.sub('[^a-zA-Z0-9./: ]','',self.df[col][i]).upper()
                    for word in bad_words:
                        if word in self.df[col][i].split(" ") or word in self.df[col][i].split(".") or word in self.df[col][i].split("/"):
                            self.df[col][i] = re.sub('[^0-9]','', self.df[col][i])
                    for word in not_too_bad_words:
                        if word in self.df[col][i].split(" "):
                            self.df[col][i] =  ' '.join(map(str, [a for a in self.df[col][i].split(" ") if a != word]))
                        if word in self.df[col][i].split(":"):
                            self.df[col][i] =  ' '.join(map(str, [a for a in self.df[col][i].split(":") 
                                                                  if a != word]))
        if (df_number == 1):
            self.df1 = self.df
            df_number = df_number + 1
        else: 
            self.df2 = self.df   
            
    def concatena_unici(self):
        a = self.df1.groupby("cleaned_species", as_index=False).agg(Conteggio = ("cleaned_species", np.size))
        onlyone1 = a[a["Conteggio"]==1]

        onlyone1 = pd.merge(onlyone1, self.df1[["cleaned_species","accessid"]] , 
                 how = 'inner',
                 left_on = "cleaned_species", 
                 right_on = "cleaned_species")
        onlyone1 = onlyone1.reset_index(drop = True)

        b = self.df2.groupby("cleaned_species", as_index=False).agg(Conteggio = ("cleaned_species", np.size))
        onlyone2 = b[b["Conteggio"]==1]


        onlyone2 = pd.merge(onlyone2, self.df2[["cleaned_species","accessid"]] , 
                 how = 'inner',
                 left_on = "cleaned_species", 
                 right_on = "cleaned_species")

        onlyone2 = onlyone2.reset_index(drop = True)

        self.single_seq = pd.merge(onlyone1,onlyone2, 
                                 how = 'inner',
                                 left_on = 'cleaned_species', 
                                 right_on = 'cleaned_species')
        
        self.df1  = self.df1[~self.df1["cleaned_species"].isin(self.single_seq['cleaned_species'].tolist())].copy()
        self.df1  = self.df1.reset_index(drop = True)
        
        self.df2 = self.df2[~self.df2["cleaned_species"].isin(single_seq['cleaned_species'].tolist())].copy()
        self.df2 = self.df2.reset_index(drop = True)  
        
    def read_concat_info(self, concat_type, organism):
        dir = os.path.join(path,"Concatenate")
        if not os.path.exists(dir):
            os.mkdir(dir)
        if platform.system() == "Darwin" or platform.system() == "Linux":
            onlyfiles = [f for f in listdir(path + "Concatenate//") if isfile(join(path + "Concatenate//", f))]
        if platform.system() == "Windows":
            onlyfiles = [f for f in listdir(path + "Concatenate\\") if isfile(join(path + "Concatenate\\", f))]
        lista_aln = [x for x in onlyfiles if x.endswith('.aln')] 
        lista_aln.sort()
        self.lista_aln = lista_aln
        lista_excel = [x for x in onlyfiles if x.endswith('.xlsx')]
        lista_excel.sort()
        self.lista_excel = lista_excel
        lista_nomi = []
        for i in range(len(lista_aln)):
            lista_nomi = lista_nomi + [lista_aln[i].split("_")[1]]
        lista_nomi.sort()
        self.lista_nomi = lista_nomi
        if concat_type == "ident":
            pass
        else:
            listadf = []
            for nomi,file_aln,file_xlsx in zip(self.lista_nomi,self.lista_aln, self.lista_excel):
                self.read_aln(alignment_file = path + "Concatenate//" + file_aln)
                self.Open_Dataframe(concat = "Y", filename = file_xlsx)
                self.merge_df_to_aln(df_type='MyDf')
                self.correcting_info(organism, df_number =1)
                listadf.append([nomi,file_aln, file_xlsx ,self.df1])
            self.listadf = listadf
    
    def add_gaps(self, Final, min_gene, min_seq,concat_type):
        print("creating a column with the summed gaps")
        lista = [x for x in list(Final.columns) if x.startswith('gap')] 
        for colonna in lista:
            Final[colonna] = np.where(Final[colonna]=="", 0, Final[colonna])
        for i in range(len(lista)-1):
            if i == 0:
                Final["Gap"] = Final[lista[i]].fillna(0)
            Final["Gap"] = Final["Gap"].fillna(0) + Final[lista[i+1]].fillna(0)
        if concat_type == "ident":
            Final['Rank'] = Final.groupby(by = ["identificativi","Species"])["Gap"].transform(
                lambda x: x.rank(method = 'first', ascending=True))
            Final = Final[Final['Rank']==1].reset_index(drop=True)
        id_col = [x for x in list(Final.columns) if x.startswith('accessid')]
        Final2 = Final.drop_duplicates(subset=id_col).copy()
        Final2 = Final2.reset_index(drop=True)
        print("checking the number of possible concatenate genes for each sample")
        Final2["Rank3"]= np.where(Final2[id_col[0]].notna(),1,0)  
        for i in range(len(id_col)-1):
            Final2["Rank3"] = Final2["Rank3"] + np.where(Final2[id_col[i+1]].notna(),1,0)
        Final2 = Final2.fillna("")
        print("Removing duplicate accessids due to multiple sequences of a single sample..")
        for col in id_col:
            Final2['Rank4'] = Final2.groupby(by = col)["Rank3"].transform(lambda x: x.rank(method = 'dense', 
                                                                                           ascending=False))
            Final2 = Final2[(Final2["Rank4"]==1) | (Final2[col]=="")]
            Final2 = Final2.reset_index(drop=True)
            Final2['Rank5'] = Final2.groupby(by = col)["Gap"].transform(lambda x: x.rank(method = 'first', 
                                                                                         ascending=True))
            Final2 =  Final2[(Final2["Rank5"]==1) | (Final2[col]=="")] 
            Final2 = Final2.reset_index(drop=True)  

        print("Keeping only concatenate with minimum", min_gene, "genes")
        Final2 = Final2[Final2["Rank3"]>=min_gene]
        Final2.reset_index(drop=True)
        
        if len(Final2)==0:
            print("no possible concatenating genes found")
            self.Concat = Final2.copy()
        else:
            if concat_type == "ident":
                print("Checking columns values:")
                print("Chosed miniumum number of sequences for each gene is: ", min_seq)
                min_null = 1 - min_seq
                numero_righe = len(Final2)
                soglia = numero_righe*min_null
                for col in id_col:
                    numero_null = len(Final2[Final2[col] == ""])
                    if numero_null > soglia:
                        gene = col.split("_")[1]
                        da_cancellare =  [x for x in list(Final2.columns) if x.endswith(gene)]  
                        Final2 = Final2.drop(da_cancellare, axis=1)     
                        print("since it has too few concatenable sequences, " + gene + " gene has been deleted")
                Final2 = Final2.reset_index(drop=True)
                self.controllo = Final2
            print("adding gaps in empty sequences")
            listaseq = [x for x in list(Final2.columns) if x.startswith('aligned')] 
            for col in listaseq: 
                if Final2[col].str.len().max() == 0: 
                    sequenza = ''.join(["-"] * int(self.memo_lunghezze[col])) 
                    Final2[col] = Final2[col].replace("",sequenza)

                elif (Final2[Final2[col]!=""][col].str.len().max() != Final2[Final2[col]!=""][col].str.len().min() ): 
                    raise ValueError('You messed up with the alignment: sequences in ' + col + ' have different lengths') 
                else: 
                    sequenza = ''.join(["-"] * int(Final2[col].str.len().max()) ) 
                    Final2[col] = Final2[col].replace("",sequenza)

            Final2["concatseq"] = Final2.loc[:,listaseq].sum(axis=1)
            Final2 = Final2.reset_index(drop=True)
            id_col = [x for x in list(Final2.columns) if x.startswith('accessid')]
            print("Concatenating the accession numbers")
            self.Concat = Final2.copy()
            for col in id_col:
                self.Concat[id_col] = self.Concat[id_col].fillna("")
            self.Concat["accessid"] = self.Concat[id_col].agg('_'.join, axis=1)
            self.Concat["accessid"] = self.Concat["accessid"].str.replace('_+', '_')
            self.Concat["accessid"] = self.Concat["accessid"].str.replace('_$', '')
            self.Concat["accessid"] = self.Concat["accessid"].str.replace('^_', '')
            
        return id_col
        
    def concat_ident(self, organism, min_gene, min_seq, concat_type):
        self.read_concat_info(concat_type = concat_type, organism = organism)
        listadf = []
        for nomi,file_aln,file_xlsx in zip(self.lista_nomi,self.lista_aln, self.lista_excel):
            self.read_aln(alignment_file = path + "Concatenate//" + file_aln)
            self.Open_Dataframe(concat = "Y", filename = file_xlsx)
            self.merge_df_to_aln(df_type='MyDf')
            self.correcting_info(organism, df_number =1)
            df1_identificativi = pd.melt(self.df1, id_vars=["accessid", "species","cleaned_species","type_material",
                                                            "country", "lat_lon", "altitude","isolation_source", 
                                                            "host", "db_xref", "collected_by", "identified_by",
                                                            "collection_date", "PCR_primers","note","sequence", 
                                                            "aligned_sequence"], 
                                               value_vars=["strain","isolate",
                                                           "specimen_voucher","clone","culture_collection"],
                                               var_name='ident_type', value_name='identificativi')
            df1_identificativi = df1_identificativi[df1_identificativi["identificativi"].notna()]
            df1_identificativi = df1_identificativi[df1_identificativi["identificativi"] != ""]
            df1_identificativi = df1_identificativi.reset_index(drop = True)  
            puntigap = []
            for i in range(len(df1_identificativi)):
                puntigap = puntigap + [df1_identificativi['aligned_sequence'][i].count('-')]
            df1_identificativi['gap'] = puntigap 
            print("creating the", nomi, "dataframe with the aligned sequence")
            listadf.append([nomi,file_aln, file_xlsx ,df1_identificativi])
        self.listadf = listadf
        
        print("merging the identifiers of the genes in a unique dataframe")
        for i in range(len(self.listadf)): 
            self.listadf[i][3] = self.listadf[i][3].add_suffix('_' + self.listadf[i][0])
            ident = [x for x in list(self.listadf[i][3].columns) if x.startswith('identificativi')] 
            self.listadf[i][3] = self.listadf[i][3].rename(columns={ ident[0] : "identificativi"})
            species = [x for x in list(self.listadf[i][3].columns) if x.startswith('cleaned_species')] 
            self.listadf[i][3] = self.listadf[i][3].rename(columns={ species[0] : "Species"})
        for i in range(len(self.listadf)-1):
            if i == 0:
                Final = self.listadf[i][3] 
            Final = pd.merge(Final, self.listadf[i+1][3],
                how = 'outer',
                left_on = ["identificativi","Species"],
                right_on = ["identificativi","Species"])
        self.Final0 = Final
        listaseq = [x for x in list(Final.columns) if x.startswith('aligned')] 
        lunghezze = []
        for col in listaseq: 
            lunghezze.append(Final[col].str.len().max())
            self.memo_lunghezze = {k:v for k,v in zip(listaseq,lunghezze)}
        id_col = self.add_gaps(Final = Final, min_gene = min_gene, min_seq = min_seq, concat_type = "ident")
        self.Concat2 = self.Concat.copy()
        self.Concat2["Rank3b"]= np.where(self.Concat2[id_col[0]] != "",1,0)  
        for i in range(1,len(id_col)):
            self.Concat2["Rank3b"] = self.Concat2["Rank3b"] + np.where(self.Concat2[id_col[i]] != "",1,0)  
        print("Keeping only concatenate sequences with minimum", min_gene, "genes")
        self.Concat2 = self.Concat2[self.Concat2["Rank3b"]>=min_gene]
        self.ConcatDf = self.Concat2.reset_index(drop = True)

        if len(self.ConcatDf) < len(self.Final0):
            print("since you have lost ", len(self.Final0) - len(self.ConcatDf),
                  " sequences, you are probably going to find some only gaps columns in Concat_Df" )
        print("ConcatDf is ready")
        self.specie_trovate_ident = list(set(list(self.ConcatDf['Species'])))
        useful_gene = [x for x in list(self.ConcatDf.columns) if x.startswith('host')]
        self.remaining_gene = [str(i).replace( 'host_', '') for i in useful_gene]
        self.Concat_df_ident = self.ConcatDf
        print(len(self.Concat_df_ident), " sequences concatenated by identificative")
        print("Concat_df_ident is ready")
  
    def concat_country(self,organism,min_gene,min_seq,concat_type):
        self.read_concat_info(concat_type="country" , organism  = organism)
        new_listadf = []
        for i in range(len(self.listadf)):
            if self.listadf[i][0] in self.remaining_gene:
                new_listadf = new_listadf + [self.listadf[i]]
        self.listadf = new_listadf.copy()
        for i in range(len(self.listadf)):
            self.listadf[i][3]["state"] = self.listadf[i][3]['country'].apply(first_word)
            self.listadf[i][3]["state"] = self.listadf[i][3]["state"].str.lower()
            Continent = pd.read_csv(self.path + "ElencoContinenti.csv") #https://github.com/AleTatti/Barcoding-Analysis/blob/main/ElencoContinenti.csv
            Continent["Nazione"] = Continent["Nazione"].str.lower()
            self.listadf[i][3]['country'] = self.listadf[i][3]['country'].str.replace('[^a-zA-Z0-9]',' ') 
            self.listadf[i][3]["raccolta"] = self.listadf[i][3]["country"].fillna('') + " " + self.listadf[i][3]["lat_lon"].fillna('')
            self.listadf[i][3]["raccolta"] = self.listadf[i][3]["raccolta"].str.lower()
            self.listadf[i][3] = pd.merge(self.listadf[i][3], Continent, 
                     how = 'left',
                     left_on = "state", 
                     right_on = "Nazione")
            self.listadf[i][3]["raccolta"] = self.listadf[i][3]["raccolta"].fillna('') + " " + self.listadf[i][3]["Continente"].fillna('')
            self.listadf[i][3]["raccolta2"]= self.listadf[i][3]["raccolta"].apply(EliminaParole) 
            self.listadf[i][3] = self.listadf[i][3][self.listadf[i][3]["raccolta2"] != "" ]
        for i in range(len(self.listadf)):
            self.listadf[i][3]= self.listadf[i][3][~self.listadf[i][3]["cleaned_species"].isin(self.specie_trovate_ident)]
        for i in range(len(self.listadf)):
            self.listadf[i][3]= self.listadf[i][3].reset_index(drop = True)
            puntigap = []
            for j in range(len(self.listadf[i][3])):
                puntigap = puntigap + [self.listadf[i][3]['aligned_sequence'][j].count('-')]
            self.listadf[i][3]['gap'] = puntigap 
        for i in range(len(self.listadf)): 
            self.listadf[i][3] = self.listadf[i][3].add_suffix('_' + self.listadf[i][0])
            species = [x for x in list(self.listadf[i][3].columns) if x.startswith('cleaned_species')] 
            self.listadf[i][3] = self.listadf[i][3].rename(columns={ species[0] : "Species"})
        for i in range(len(self.listadf)-1):
            if i == 0:
                Final = self.listadf[i][3] 
            Final = pd.merge(Final, self.listadf[i+1][3],
                how = 'outer',
                left_on = ["Species"],
                right_on = ["Species"])       
        Final = Final.reset_index(drop = True)
        Final = Final.fillna("")
        id_raccolta2 = [x for x in list(Final.columns) if x.startswith('raccolta2')]
        lista_punteggi = []
        for index in range(len(Final)):
            punteggio = 0
            for i, colonna in zip(range(len(id_raccolta2)),id_raccolta2):
                for j, colonna2 in zip(range(i+1,len(id_raccolta2)), id_raccolta2[i+1:]):
                    punteggiop = len(set(Final.loc[index,colonna]).intersection(set(Final.loc[index,colonna2])))
                    punteggio = punteggio +  punteggiop
            lista_punteggi.append(punteggio) 
        Final["punteggio"] = lista_punteggi   
        df_country0 = Final.copy()
        df_country0['Rank_country'] =  df_country0.groupby(by = ["Species"])["punteggio"].transform(
            lambda x: x.rank(method = 'first', ascending=False))
        Final = df_country0[(df_country0['Rank_country'] == 1) & (df_country0['punteggio'] > 0) ].copy()
        print("the merged df_country is ready")
        ListaNonTrovati = df_country0[(df_country0['Rank_country'] == 1) & (df_country0['punteggio'] == 0)]["Species"].tolist()
        self.country_non_trovati = ListaNonTrovati
        self.not_merged = df_country0[df_country0["Species"].isin(ListaNonTrovati)]
        self.Final = Final.copy()
        id_col = self.add_gaps(Final = Final, min_gene = min_gene, min_seq = min_seq, concat_type = "country")
        print(id_col)
        self.Concat_df_country = self.Concat.copy() 
        print(len(self.Concat_df_country), " chimaeric species concatenated by sampling place")
        print("Concat_df_country is ready")
        
    def concat_host(self,organism,min_gene,min_seq,concat_type):
        self.read_concat_info(concat_type="host" , organism  = organism)
        new_listadf = []
        for i in range(len(self.listadf)):
            if self.listadf[i][0] in self.remaining_gene:
                new_listadf = new_listadf + [self.listadf[i]]
        self.listadf = new_listadf.copy()
        for i in range(len(self.listadf)):
            self.listadf[i][3] = self.listadf[i][3][~self.listadf[i][3]["cleaned_species"].isin(
                self.specie_trovate_ident) ]
        for i in range(len(self.listadf)):
            self.listadf[i][3]['host'] = self.listadf[i][3]['host'].str.replace('[^a-zA-Z0-9]',' ')
            
        for i in range(len(self.listadf)):
            self.listadf[i][3] = self.listadf[i][3].reset_index(drop = True)
            puntigap = [] 
            for j in range(len(self.listadf[i][3])):
                puntigap = puntigap + [self.listadf[i][3]['aligned_sequence'][j].count('-')]
            self.listadf[i][3]['gap'] = puntigap 
        for i in range(len(self.listadf)): 
            self.listadf[i][3] = self.listadf[i][3].add_suffix('_' + self.listadf[i][0])
            species = [x for x in list(self.listadf[i][3].columns) if x.startswith('cleaned_species')] 
            self.listadf[i][3] = self.listadf[i][3].rename(columns={ species[0] : "Species"})
            host = [x for x in list(self.listadf[i][3].columns) if x.startswith('host')] 
            self.listadf[i][3] = self.listadf[i][3].rename(columns={ host[0] : "host"})

        for i in range(len(self.listadf)-1):
                if i == 0:
                    Final = self.listadf[i][3] 
                Final = pd.merge(Final, self.listadf[i+1][3],
                    how = 'outer',
                    left_on = ["host", "Species"],
                    right_on = ["host", "Species"]) 
        
        if concat_type == "best_chimaeras": 
            Final = Final[Final["Species"].isin(self.country_non_trovati)]
            Final = Final.reset_index(drop = True)
            self.controllo_Fin = Final

        id_col = self.add_gaps(Final = Final, min_gene = min_gene, min_seq = min_seq, concat_type = "host")
        self.Concat_df_host = self.Concat.copy()
        print(len(self.Concat_df_host), " chimaeric species concatenated by host")
        print("Concat_df_host is ready")
        
        only_chimaeras = Final[~Final["Species"].isin(self.Concat_df_host['Species'])]
        
        '''[~self.Concat_df_host["Species"].isin(self.country_non_trovati)].copy()'''
        self.only_chimaeras = only_chimaeras.reset_index(drop = True)
        print(len(list(set(list(self.only_chimaeras["Species"])))), " species will be randomly concatenated")
        
               
    def concat_random_chimaeras(self,organism,min_gene, min_seq,concat_type):
        Final =  self.only_chimaeras
        Final['Rank'] = Final.groupby(by = ["Species"])["Gap"].transform(
            lambda x: x.rank(method = 'first', ascending=True))
        Final = Final[Final["Rank"]==1].copy()
        Final = Final.reset_index(drop = True)   
        id_col = self.add_gaps(Final = Final, min_gene = min_gene, min_seq = min_seq, concat_type = "best_chimaeras")
        Final = self.Concat.copy()
        self.Concat_df_chimaeras = Final.copy()
        print(len(self.Concat_df_chimaeras), " chimaeric species randomly concatenated")
        print("Concat_df_chimaeras is ready")
    
    def concat_sequences(self,organism,min_gene,min_seq,concat_type):
        if concat_type == "ident":
            self.concat_ident(organism,min_gene,min_seq,concat_type="ident")
            self.df_Concat = self.Concat_df_ident
        else:
            if concat_type == "best_chimaeras":
                self.concat_ident(organism,min_gene,min_seq,concat_type="ident")
                self.concat_country(organism,min_gene,min_seq,concat_type="country")
                self.concat_host(organism,min_gene,min_seq,concat_type="best_chimaeras")
                self.concat_random_chimaeras(organism,min_gene,min_seq,concat_type="best_chimaeras")
                df_Concat = pd.concat([self.Concat_df_ident,
                                       self.Concat_df_country, 
                                       self.Concat_df_host, 
                                       self.Concat_df_chimaeras],
                      ignore_index = True)

            if concat_type == "host":
                self.concat_ident(organism,min_gene,concat_type="ident", min_seq=0.5)
                self.concat_host(organism,min_gene,concat_type="host", min_seq=0.5)
                df_Concat = pd.concat([self.Concat_df_ident, 
                                       self.Concat_df_host], 
                      ignore_index = True)
            if concat_type != "host" and concat_type != "country" and concat_type != "best_chimaeras":
                print("concat_type could be: ident, host, country, and best_chimaeras")

            df_Concat = df_Concat.reset_index(drop=True)
            self.df_Concat = df_Concat.copy()
        print("a total of ", len(self.df_Concat), " sequences have been concatenated")
        print("df_Concat is ready")
        self.Save_fasta(cleaning_level=1,n=0.9,G=0.9,Dataframe=self.df_Concat,out="Concat_df",header="default")
            
    def Help(self, key):
        if key == "read_file":
            prin('''read_file(filename,tipologia) where tipologia could be: "list" (meaning a list separated by \\n) 
            or "fasta" meaning fasta file''')
        
        elif key == "filename":
            print("filename could be:") 
            print("a txt file with a list separated by \\n")
            print("a fasta file composed by >header \n sequence")
            print("a list separated by a comma without spaces")
        
        elif key == "tipologia":
                  print("tipologia could be:")
                  print('''"list" (meaning a list separated by \\n)''')
                  print('''"fasta (meaning a fasta file composed by >header \\n sequence)"''')
        
        elif key == "transform_file_in_filter":
            print("transform_file_in_filter(filename, gene, rank_type) where:")
            print('''"gene" is a string with the name of the marker you are going to search, 
            i.e. "Internal Transcribed Spacer"''')
            print('''"rank_type" has the default 1 meaning searching fot the genus or superior rank, 
            otherwise it is going to search for each species''')

        
        elif key == "gene":
                  print('''"gene" is a string with the name of the marker you are going to search, i.e. "Internal Transcribed Spacer"''')
                  print('''Be careful: "Internal Transcribed Spacer" is different from "ITS"
                  you need to use the full name in the keyword search but the initials for the search inside a genome''')
        
        elif key == "rank_type": 
                  print('''"rank_type" mean if you need to search for the species or higher taxa''')
                  print("It has the default 1 meaning searching fot the genus or superior rank such as family")
                  print("If youMyDf put something else it is going to search for the species")
        
        
        elif key == "find_accessid":
                  print("find_accessid(par_retmax) where par_retmax is the maximum number of sequence you want to obtain")
        
        elif key == "execute_blast":
                  print('''execute_blast(query,threshold,csv)''')
                  print('''"query" is a fasta file with the sequence you need to blast''')
                  print('''default=90: "threshold" is the identity threshold you need to keep''')
                  print('''default="none": "csv" is the file obtained from blast; 
                  \nthe "none" default means is the one automatically obtained by the run_blast function''')
        elif key == "query":
                  print('''"query" is a fasta file with the sequence you need to blast''')
        elif key == "threshold":
                  print('''default=90: "threshold" is the identity threshold you need to keep''')
        elif key == "csv":
                  print('''default="none": "csv" is the file obtained automatically from blast''')
                  print('''"none" means the one obtained authomatically by the run_blast function''')
                  print('''you can use a csv directely obtained from NCBI blast''')
        elif key == "find_from_accessid":
            print('''find_from_accessid(query, threshold, csv, filename, id_type, gene, rank_type)''')
        elif key == "id_type":
            print('''"id_type" is the input id; it could be:''')
            print('''"total": meaning that the search is going to use both keyword and blast search to create a id list''')
            print('''"B": meaning that the search is going to use blast to create a id list''')
            print('''"K": meaning that the search is going to use the keyword search to create a id list''')
            print('''"file": meaning that the search uses a input file (filename) with a list of ids separated by "\\n" ''')
            print('''"list": meaning that the search uses a list separated by comma given in filename''')
        elif key == "df_type":
            print('''"df_type" could be: \n "df" meaning the authomatic one obtained from the search \n "MyDf" meaning the one uploaded''')
        elif key == "concat_ident":
            print('''concat_ident(organism,above,min_seq=0.5) \n
            is the concatenation function for identification. \n
            it requires the alignment file and the dataframe with the info for at least two genes. \n
            They must be inside the directory called "Concatenate" and they must be named in the following way: \n
            organism_GeneName.xlsx, organism_GeneName.aln \n
            ie. Leucoagaricus_ITS.xlsx, Leucoagaricus_ITS.aln ; Leucoagaricus_LSU.xlsx, Leucoagaricus_LSU.aln \n
            ''')
        elif key == "organism":
            print('''"organism" can be "Fungi" or "Wolbachia" or "other" needed because of the peculiar way the identification info are written in those groups''')
        elif key == "min_gene":
            print('''min_gene is the minimum number of genes I want to concatenate \n 
                  i.e. 3 means that I'm discarding all the speciment with only 2 concatenate genes''')
    
timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
##creo la classe BarcodingGap che eredita i metodi di Create_Dataset
class BarcodingGap(Create_Dataset):
       
    def __init__(self, path, blastn_exe, mafft_exe,iqtree): 
        super().__init__(path, blastn_exe, mafft_exe,iqtree)
        self.path = path
        self.blastn_exe = blastn_exe
        self.mafft_exe = mafft_exe

    def read_file_records(self,filename):
        self.filename = filename
        self.records = SeqIO.parse(self.path + self.filename, "fasta")
        
    def convert_file(self): 
        self.file_phylip = self.filename + ".phylip"
        self.count = SeqIO.write(self.records, self.path + self.file_phylip, "phylip")
        
    def create_matrix(self): 
        aln = AlignIO.read(open(self.path + self.file_phylip), 'phylip') 
        calculator = DistanceCalculator() 
        self.dm = calculator.get_distance(aln) 
        
    def create_dataframe(self):
        lista_nomi = self.dm.names
        matrice = self.dm.matrix
        MyList = []
        for i in range(len(lista_nomi)):
            for j in range(len(lista_nomi)):
                if j>i:     
                    MyList.append([lista_nomi[i],i,lista_nomi[j],j])
        MyDf = pd.DataFrame(MyList, columns=['Sequenza1','Posizione1', 'Sequenza2','Posizione2'])            
        Distanza = []
        for i in range(len(MyDf)):
            Distanza.append(matrice[MyDf.iloc[i,3]][MyDf.iloc[i,1]])
        MyDf2 = pd.DataFrame(Distanza, columns=['Distanza'])
        MyDf = pd.concat([MyDf, MyDf2],
          axis=1)
        MyDf = MyDf.drop(['Posizione1'],axis=1)
        MyDf = MyDf.drop(['Posizione2'],axis=1) 
        self.MyDf = MyDf    
        
    def extract_info_from_fasta(self,filename):
        self.filename = filename
        o = open(self.path+self.filename, "r")
        alignment = o.read()
        o.close

        df = pd.DataFrame(columns = ['accessid','species','other_info','aligned_sequence'] )
        dati_splittati = alignment.split(">")
        dati_splittati.pop(0)

        data  = []
        for i in range(len(dati_splittati)):

            header = dati_splittati[i].split("\n")[0]
            prefix = header.split("_")[0]
            seq = dati_splittati[i].split("\n")[1]
            
            if filename.startswith("Concat"):
                if(len(prefix)==0):
                    data = data + [{'accessid': header.split("_")[1] + "_" + header.split("_")[2],
                                                'species': header.split("_")[-2] + "_" + header.split("_")[-1],
                                                'other_info': "",
                                       'aligned_sequence': seq}]
                else:
                    data = data + [{'accessid': header.split("_")[0],
                                                'species': header.split("_")[-2] + "_" + header.split("_")[-1],
                                                'other_info': "",
                                       'aligned_sequence': seq}]
            else:
                if(len(prefix) == 0): 
                    if(len(header.split("_")[2]) < 3):  
                        if((header.split("_")[4] == 'uncultured') or (header.split("_")[4] == "Uncultured") 
                                                                  or (header.split("_")[4] == "UNVERIFIED")): 
                            if(header.split("_")[6] == 'cf.' or header.split("_")[6] == 'aff.' or header.split("_")[6] == 'nr.'): 
                                if len(header.split("_"))>7: 
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                'species': header.split("_")[5]+ "_" + header.split("_")[7],
                                                'other_info': '_'.join(header.split("_")[8:]),
                                       'aligned_sequence': seq}]
                                else:    
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                    'species': header.split("_")[5]+ "_" + header.split("_")[7],
                                           'aligned_sequence': seq}]
                            else:
                                if len(header.split("_"))>6:
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                    'species': header.split("_")[5]+ "_" + header.split("_")[6],
                                                    'other_info': '_'.join(header.split("_")[7:]),
                                           'aligned_sequence': seq}]  
                                else:
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                'species': header.split("_")[5]+ "_" + header.split("_")[6],
                                       'aligned_sequence': seq}]  

                            if(header.split("_")[6][0].isupper()): 
                                data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                'species': header.split("_")[5] + "_sp.",
                                                'other_info': '_'.join(header.split("_")[6:]),
                                       'aligned_sequence': seq}]    
                        else:    
                            if(header.split("_")[5] == 'cf.' or header.split("_")[5] == 'aff.'): 
                                if len(header.split("_"))>6:
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                    'species': header.split("_")[4]+ "_" + header.split("_")[6],
                                                    'other_info': '_'.join(header.split("_")[7:]),
                                           'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                    'species': header.split("_")[4]+ "_" + header.split("_")[6],
                                           'aligned_sequence': seq}]
                            else:
                                if(header.split("_")[5][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                'species': header.split("_")[4] + "_sp.",
                                                'other_info': '_'.join(header.split("_")[5:]),
                                       'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[2] + "_" + header.split("_")[3],
                                                    'species': header.split("_")[4]+ "_" + header.split("_")[5],
                                           'aligned_sequence': seq}]
                    else:
                        if((header.split("_")[3] == 'uncultured') or (header.split("_")[3] == "Uncultured" 
                                                                      or (header.split("_")[3] == "UNVERIFIED"))): 
                            if(header.split("_")[5] == 'cf.' or header.split("_")[5] == 'aff.'):
                                if len(header.split("_"))>6:
                                    data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[4]+ "_" + header.split("_")[6],
                                                    'other_info': '_'.join(header.split("_")[7:]),
                                                    'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[4]+ "_" + header.split("_")[6],
                                                    'aligned_sequence': seq}]
                            else:
                                if (header.split("_")[5][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[2], 
                                                'species':header.split("_")[4] + "_sp.",
                                                'other_info': '_'.join(header.split("_")[5:]),
                                                'aligned_sequence': seq}]
                                else:    
                                    if len(header.split("_"))>5:
                                        data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[4]+ "_" + header.split("_")[5],
                                                    'other_info': '_'.join(header.split("_")[6:]),
                                                    'aligned_sequence': seq}]
                                    else:
                                        data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[4]+ "_" + header.split("_")[5],
                                                    'aligned_sequence': seq}]
                        else:
                            if(header.split("_")[4] == 'cf.' or header.split("_")[4] == 'aff.'):
                                if len(header.split("_"))>5:
                                    data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[3]+ "_" + header.split("_")[5],
                                                    'other_info': '_'.join(header.split("_")[6:]),
                                                    'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[2], 
                                                    'species':header.split("_")[3]+ "_" + header.split("_")[5],
                                                    'aligned_sequence': seq}]
                            else: 
                                if(header.split("_")[4][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[2], 
                                                'species':header.split("_")[3] + "_sp.",
                                                'other_info': '_'.join(header.split("_")[4:]),
                                                'aligned_sequence': seq}]
                                else:
                                    if len(header.split("_"))>4:
                                        data = data + [{'accessid': header.split("_")[2], 
                                                        'species':header.split("_")[3]+ "_" + header.split("_")[4],
                                                        'other_info': '_'.join(header.split("_")[5:]),
                                                        'aligned_sequence': seq}]
                                    else:
                                        data = data + [{'accessid': header.split("_")[2], 
                                                        'species':header.split("_")[3]+ "_" + header.split("_")[4],
                                                        'aligned_sequence': seq}]
                else: 
                    if(len(prefix) < 3):
                        if((header.split("_")[3] == 'uncultured') or (header.split("_")[3] == "Uncultured")
                                                                   or (header.split("_")[3] == "UNVERIFIED")): 
                            if(header.split("_")[5] == 'cf.' or header.split("_")[5] == 'aff.'):
                                if len(header.split("_"))>5:
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                    'species': header.split("_")[3]+ "_" + header.split("_")[5],
                                                    'other_info': '_'.join(header.split("_")[6:]),
                                                    'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                'species': header.split("_")[3]+ "_" + header.split("_")[5],
                                         'aligned_sequence': seq}]
                            else:
                                if(header.split("_")[4][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                    'species': header.split("_")[3] + "_sp.",
                                                    'other_info': '_'.join(header.split("_")[4:]),      
                                                    'aligned_sequence': seq}] 
                                else:
                                    if len(header.split("_"))>4:
                                        data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                        'species': header.split("_")[3]+ "_" + header.split("_")[4],
                                                        'other_info': '_'.join(header.split("_")[5:]),   
                                                        'aligned_sequence': seq}] 
                                    else:
                                        data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                        'species': header.split("_")[3]+ "_" + header.split("_")[4],
                                                        'aligned_sequence': seq}] 
                        else:
                            if(header.split("_")[3] == 'cf.' or header.split("_")[3] == 'aff.'):
                                if len(header.split("_"))>4:
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                'species': header.split("_")[2]+ "_" + header.split("_")[4],
                                                'other_info': '_'.join(header.split("_")[5:]),  
                                                'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                    'species': header.split("_")[2]+ "_" + header.split("_")[4],
                                                    'aligned_sequence': seq}]
                            else:
                                if(header.split("_")[3][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                    'species': header.split("_")[2] + "_sp.",
                                                    'other_info': '_'.join(header.split("_")[3:]),  
                                                    'aligned_sequence': seq}]
                                else:    
                                    if len(header.split("_"))>3:
                                        data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                        'species': header.split("_")[2]+ "_" + header.split("_")[3],
                                                        'other_info': '_'.join(header.split("_")[4:]),  
                                                        'aligned_sequence': seq}]
                                    else:
                                        data = data + [{'accessid': header.split("_")[0] + "_" + header.split("_")[1],
                                                        'species': header.split("_")[2]+ "_" + header.split("_")[3],
                                                        'aligned_sequence': seq}]
                    else: 
                        if((header.split("_")[1] == 'uncultured') or (header.split("_")[1] == "Uncultured")
                                                                   or (header.split("_")[1] == "UNVERIFIED")): 
                            if (len(header.split("_")) >3): 
                                if(header.split("_")[3] == 'cf.' or header.split("_")[3] == 'aff.'):
                                    if len(header.split("_"))>4:
                                        data = data + [{'accessid': header.split("_")[0],
                                                        'species': header.split("_")[2]+ "_" + header.split("_")[4],
                                                        'other_info': '_'.join(header.split("_")[5:]),  
                                                        'aligned_sequence': seq}]
                                    else:
                                        data = data + [{'accessid': header.split("_")[0],
                                                        'species': header.split("_")[2]+ "_" + header.split("_")[4],
                                                        'aligned_sequence': seq}]
                                else:
                                    if(header.split("_")[3][0].isupper()): 
                                        data = data + [{'accessid': header.split("_")[0],
                                                        'species': header.split("_")[2] + "_sp.",
                                                        'other_info': '_'.join(header.split("_")[3:]),  
                                                        'aligned_sequence': seq}]
                                    else:
                                        if len(header.split("_"))>3:
                                            data = data + [{'accessid': header.split("_")[0],
                                                            'species': header.split("_")[2]+ "_" + header.split("_")[3],
                                                            'other_info': '_'.join(header.split("_")[4:]),  
                                                            'aligned_sequence': seq}]
                                        else:
                                            data = data + [{'accessid': header.split("_")[0],
                                                            'species': header.split("_")[2]+ "_" + header.split("_")[3],
                                                            'aligned_sequence': seq}]
                            if (len(header.split("_")) ==3): 
                                    data = data + [{'accessid': header.split("_")[0],
                                                    'species': header.split("_")[2]+ "_sp.",
                                                    'aligned_sequence': seq}]       
                        else:
                            if(header.split("_")[1] == 'cf.' or header.split("_")[2] == 'aff.'):
                                if len(header.split("_"))>3:
                                    data = data + [{'accessid': header.split("_")[0],
                                                    'species': header.split("_")[1]+ "_" + header.split("_")[3],
                                                    'other_info': '_'.join(header.split("_")[4:]), 
                                                    'aligned_sequence': seq}]
                                else:
                                    data = data + [{'accessid': header.split("_")[0],
                                                    'species': header.split("_")[1]+ "_" + header.split("_")[3],
                                                    'aligned_sequence': seq}]
                            else:
                                if(header.split("_")[2][0].isupper()):
                                    data = data + [{'accessid': header.split("_")[0],
                                                    'species': header.split("_")[1] + "_sp.",
                                                    'other_info': '_'.join(header.split("_")[2:]), 
                                                    'aligned_sequence': seq}]
                                else:
                                    if len(header.split("_"))>2:
                                        data = data + [{'accessid': header.split("_")[0],
                                                        'species': header.split("_")[1]+ "_" + header.split("_")[2],
                                                        'other_info': '_'.join(header.split("_")[3:]), 
                                                        'aligned_sequence': seq}]
                                    else:
                                        data = data + [{'accessid': header.split("_")[0],
                                                        'species': header.split("_")[1]+ "_" + header.split("_")[2],
                                                        'aligned_sequence': seq}]
        df = df.append(data,ignore_index=True)
        self.df = df   
        
    def merge_df_mydf(self):
        df2 = self.df[self.df["accessid"].str.len()>10].copy()
        df2 = df2.reset_index(drop=True)
        df2["accessid_primidieci"] = df2.accessid.str[0:10]
        self.MyDf["sequenza1_primiotto"] = self.MyDf.Sequenza1.str[0:8]
        self.MyDf["sequenza2_primiotto"] = self.MyDf.Sequenza2.str[0:8]
        df3 = self.df[self.df["accessid"].str.len()<10].copy()
        df3 = df3.reset_index(drop=True)
        df3["accessid_primiotto"] = df3.accessid.str[0:8]
        merge0 = pd.merge(self.MyDf, df2, how='left', left_on="Sequenza1", right_on="accessid_primidieci")
        self.MyDf['Sequenza1'] = np.where(merge0["accessid"].notna(),
                                            merge0["accessid"],
                                            self.MyDf['Sequenza1'])
        merge0 = pd.merge(self.MyDf, df2, how='left', left_on="Sequenza2", right_on="accessid_primidieci")
        self.MyDf['Sequenza2'] = np.where(merge0["accessid"].notna(),
                                          merge0["accessid"],
                                          self.MyDf['Sequenza2'])
        merge0 = pd.merge(self.MyDf, df3, how='left', left_on="sequenza1_primiotto", right_on="accessid_primiotto")
        self.MyDf['Sequenza1'] = np.where(merge0["accessid"].notna() ,
                                            merge0["accessid"],
                                            self.MyDf['Sequenza1'])
        merge0 = pd.merge(self.MyDf, df3, how='left', left_on="sequenza2_primiotto", right_on="accessid_primiotto")
        self.MyDf['Sequenza2'] = np.where(merge0["accessid"].notna(),
                                          merge0["accessid"],
                                          self.MyDf['Sequenza2']) 
        self.MyDf = self.MyDf.drop(['sequenza1_primiotto'],axis=1)
        self.MyDf = self.MyDf.drop(['sequenza2_primiotto'],axis=1)

        merge1 = pd.merge(self.MyDf, self.df, how='inner', left_on="Sequenza1", right_on="accessid")
        self.merge1 = merge1
        merge1 = merge1[["Sequenza1","Sequenza2","Distanza","species","other_info"]]
        self.merge1_0 = merge1
        merge1 = merge1.rename(columns={"species": "Specie1", "other_info":"other_info1"})
        self.merge1_1 = merge1
        merge2  = pd.merge(merge1, self.df, how='inner', left_on="Sequenza2", right_on="accessid")
        merge2 = merge2[["Sequenza1","Sequenza2","Distanza","Specie1","species","other_info1","other_info"]]
        merge2 = merge2.rename(columns={"species": "Specie2", "other_info":"other_info2"})
        self.merge1_2 = merge2
        if(len(self.MyDf) != len(merge2)):
            print("ALERT Check Access ID or you are going to lose data")
        self.merge2 = merge2
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        self.merge2.to_csv("merge2_dataframe_" + timestamp + ".csv")

    def save_excel(self,dataframe): 
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        dataframe.to_excel(self.path + "distance_matrix_" + timestamp + ".xlsx", index = False)
        
    def create_dataframe_matrix(self,filename):
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + " reading the input file")
        self.read_file_records(filename)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + " converting file fasta in phylip")
        self.convert_file()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" converted %i records" % self.count)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" creating the distance matrix ...")
        self.create_matrix()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" creating the distance dataframe")
        self.create_dataframe()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" modifying unprecise definitions such as affinis or confer")
        self.extract_info_from_fasta(filename=filename)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" correcting truncated names")
        self.merge_df_mydf()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        if (len(self.merge2) < 1048576):
            print(timestamp +" saving distance dataframe in excel")
            self.save_excel(self.merge2)
        else:
            print(timestamp +" since the dataframe is too long: saving distance dataframe in csv")
            self.merge2.to_csv(timestamp + "_merge2.csv")
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp +" done")
        
    def process_matrix(self, rank, kingdom, interesting_species):
        merge2 = self.merge2.copy()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        if rank == "Population":
            
            print(timestamp + ' reading the interesting populations input')
            interest_specific_name = interesting_species.split(",")
            merge2["DistanzaStringa"] = merge2["Distanza"]
            merge2["DistanzaStringa"] = merge2["DistanzaStringa"].astype(str)
            merge2 = merge2.rename(columns={"other_info1": "id1", "other_info2":"id2"})  
            self.controllopop0 = merge2
            merge2.replace({"id1": r'_'}, {'id1': ' '}, regex=True, inplace=True)
            self.controllopop1 = merge2
            merge2.replace({"id2": r'_'}, {'id2': ' '}, regex=True, inplace=True)
            self.controllopop2 = merge2
            merge2.reset_index(drop=True)
                
        if rank == "Species":
            print(timestamp + ' reading the interesting species input')
            interest_genus = []
            interest_specific_name = []
            interesting_species = interesting_species.split(",")
            for i in range(len(interesting_species)):
                interest_genus = interest_genus + [interesting_species[i].split(" ")[0]]
                interest_genus = list(set(interest_genus))
                interest_specific_name = interest_specific_name + [interesting_species[i].split(" ")[1]]
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            print(timestamp + ' assigning the chosen rank')
            merge2[["genus1","id1"]] = merge2.Specie1.str.split("_",expand=True)
            merge2[["genus2","id2"]] = merge2.Specie2.str.split("_",expand=True)
            merge2 = merge2[ (merge2["genus1"] != "Unidentified") & (merge2["genus2"] != "Unidentified") ]
            merge2.reset_index(drop=True)
            merge2 = merge2[ (merge2["id1"] != "sp.") & (merge2["id2"] != "sp.") ] 
            merge2.reset_index(drop=True)
            self.merge2_0 = merge2
            if kingdom == "animals": 
                for x in ['oidea','idae','inae']:
                    merge2 = merge2[~merge2["genus1"].str.endswith(x)]
                    merge2.reset_index(drop=True)
                    merge2 = merge2[~merge2["genus2"].str.endswith(x)]
                    merge2.reset_index(drop=True)
            if kingdom == "fungi":
                for x in ['mycota','mycotina','mycetes', 'mycetidae', 'ales', 'ineae', 'aceae']:
                    merge2 = merge2[~merge2["genus1"].str.endswith(x)]
                    merge2 = merge2[~merge2["genus2"].str.endswith(x)]
                    merge2.reset_index(drop=True)
                merge2 = merge2[ (merge2["genus1"] != "uncultured") & (merge2["genus2"] != "uncultured") ]
                merge2.reset_index(drop=True)
                merge2 = merge2[ (merge2["genus1"] != "fungal") & (merge2["genus2"] != "fungal") ]
                merge2.reset_index(drop=True)
            if kingdom == "plants":
                for x in ['spermae','phyta','phytina','phyceae','phycidae', 'anae', 'ales', 'ineae', 'aria','acea','aceae', 'oideae']:
                    merge2 = merge2[~merge2["genus1"].str.endswith(x)]
                    merge2 = merge2[~merge2["genus2"].str.endswith(x)]
                    merge2.reset_index(drop=True)     
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')    
            print(timestamp + ' inverting column values too keep consistency in the comparison')   
            merge2["DistanzaStringa"] = merge2["Distanza"]
            merge2["DistanzaStringa"] = merge2["DistanzaStringa"].astype(str)
            merge2.reset_index(drop=True)
        merge2['Sequenza1'],merge2['id1'],merge2['Specie1'],merge2['Sequenza2'],merge2['id2'],merge2['Specie2'] = np.where(merge2['id1']>merge2['id2'],
             (merge2['Sequenza2'],merge2['id2'],merge2['Specie2'],merge2['Sequenza1'],merge2['id1'],merge2['Specie1']),
             (merge2['Sequenza1'],merge2['id1'],merge2['Specie1'],merge2['Sequenza2'],merge2['id2'],merge2['Specie2']))
        merge2.reset_index(drop=True)   
            
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + ' calculating the intraspecific and interspecific ranges')
        merge2["id1_id2"] = merge2["id1"] + "_" + merge2["id2"]
        self.merge2_1 = merge2
        intra = merge2[ (merge2["id1"] == merge2["id2"]) & (~merge2["id1"].isin(interest_specific_name)) ]["Distanza"].tolist()
        self.intra = intra
        verify_intra = merge2[ (merge2["id1"] == merge2["id2"]) & (~merge2["id1"].isin(interest_specific_name)) ]
        self.verity_intra = verify_intra.reset_index(drop=True)
        inter = merge2[ (merge2["id1"] != merge2["id2"])&(~merge2["id1"].isin(interest_specific_name))&
                       (~merge2["id2"].isin(interest_specific_name))]["Distanza"].tolist()
        self.inter = inter
        self.impossible_inter = merge2[merge2["id1"] != merge2["id2"]]
        self.impossible_inter = self.impossible_inter[self.impossible_inter["Distanza"]==0]
        verify_inter = merge2[ (merge2["id1"] != merge2["id2"]) & 
                              (~merge2["id1"].isin(interest_specific_name)) & (~merge2["id2"].isin(interest_specific_name)) ]
        self.verify_inter = verify_inter.reset_index(drop=True)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + ' choosing the interesting combinations')
        merge3 = merge2.copy()
        merge3 = merge3[(merge2["id1"].isin(interest_specific_name)) | (merge3["id2"].isin(interest_specific_name)) ]   
        merge3.reset_index(drop=True)
        merge3 = merge3.groupby(['id1_id2'])['DistanzaStringa'].apply(';'.join).reset_index()
        merge3["DistanzaStringaLista"]=merge3["DistanzaStringa"].apply(tolist)
        merge3.reset_index(drop=True)
        merge4 = merge3.set_index('id1_id2')
        merge4 = merge4.drop(['DistanzaStringa'],axis=1)
        self.merge4 = merge4.reset_index(drop=True)
        dict = merge4.to_dict()
        dict = dict['DistanzaStringaLista']
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')     
        print(timestamp + " calculating the barcoding gap")
        upper_quartile_intra = np.percentile(intra, 75)
        lower_quartile_intra = np.percentile(intra, 25)
        lower_wisker_intra = 0.0
        upper_wisker_intra = upper_quartile_intra + 1.5 * (upper_quartile_intra - lower_quartile_intra)
        if upper_wisker_intra > max(intra):
            upper_wisker_intra = max(intra)
        upper_quartile_inter = np.percentile(inter, 75)
        lower_quartile_inter = np.percentile(inter, 25)
        lower_wisker_inter = lower_quartile_inter - 1.5 * (upper_quartile_inter - lower_quartile_inter)
        if lower_wisker_inter < min(inter):
            lower_wisker_inter = min(inter)
        upper_wisker_inter = upper_quartile_inter + 1.5 * (upper_quartile_inter - lower_quartile_inter)
        if upper_wisker_inter > max(inter):
            upper_wisker_inter = max(inter)
        if upper_wisker_intra < lower_wisker_inter:
            print("the barcoding gap is:", (upper_wisker_intra , lower_wisker_inter))
        else:
            print("the barcoding gap seems not present")
        self.merge22 = merge2
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + ' saving intraspecific and interspecific outliers')
        self.strange = merge2[(merge2["Distanza"] > upper_wisker_intra) &  (merge2["Distanza"] < lower_wisker_inter)]
        outlier_intra = merge2[(merge2["id1"] == merge2["id2"]) &  (merge2["Distanza"] > upper_wisker_intra)]
        self.outlier_intra = outlier_intra.reset_index(drop=True)
        lower_outlier_inter = merge2[(merge2["id1"] != merge2["id2"]) &  (merge2["Distanza"] < lower_wisker_inter)]
        if len(lower_outlier_inter)==0:
            lower_outlier_inter = self.impossible_inter
        self.lower_outlier_inter = lower_outlier_inter.reset_index(drop=True)
        upper_outlier_inter = merge2[(merge2["id1"] != merge2["id2"]) &  (merge2["Distanza"] > upper_wisker_inter )]
        self.upper_outlier_inter = upper_outlier_inter.reset_index(drop=True)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + ' removing outliers from the intraspecific and interspecific range')
        self.final_intra_range = merge2[(merge2["id1"] == merge2["id2"]) &  
                                        (merge2["Distanza"] < upper_wisker_intra ) & (~merge2["id1"].isin(interest_specific_name))]["Distanza"].tolist() 
        self.final_inter_range = merge2[(merge2["id1"] != merge2["id2"]) &  
                                        (merge2["Distanza"] > lower_wisker_inter ) &  (merge2["Distanza"]  < upper_wisker_inter ) & 
                                        (~merge2["id1"].isin(interest_specific_name)) & (~merge2["id2"].isin(interest_specific_name))]["Distanza"].tolist() 
        self.upper_wisker_intra = upper_wisker_intra
        self.lower_wisker_inter = lower_wisker_inter
        self.upper_wisker_inter = upper_wisker_inter   
        self.interest_specific_name = interest_specific_name 
        self.intra = intra
        self.inter = inter
        self.newmerge = merge2
        self.dict = dict
        self.merge3 = merge3
        self.verify_intra = verify_intra
        self.verify_inter = verify_inter
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp + 'done')
        
    def find_outlier(self, search_type, filename):
        lista_outliers = [self.outlier_intra, self.lower_outlier_inter, self.upper_outlier_inter]
        print("Saving the outliers in a excel file")
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        with pd.ExcelWriter(path+ 'output_outliers_' + timestamp + '.xlsx') as writer:
            self.outlier_intra.to_excel(writer, sheet_name='outlier_intra')
            self.lower_outlier_inter.to_excel(writer, sheet_name='lower_outlier_inter')
            self.upper_outlier_inter.to_excel(writer, sheet_name='upper_outlier_inter')
        print("studying the outliers..")
        for i in range(len(lista_outliers)):
            find_seq0 = pd.melt(lista_outliers[i], 
                                id_vars=[], 
                                value_vars=["Sequenza1","Sequenza2"],
                                var_name='ident_type', 
                                value_name='id_outlier')
            find_seq1 = find_seq0["id_outlier"].value_counts().to_frame(name="occurrence") 
            find_seq1['accessid'] = find_seq1.index 
            find_seq1 = find_seq1.reset_index().drop("index",1) 
            find_seq2 = find_seq1[find_seq1["occurrence"]>1].copy()
            accessid = list(find_seq2["accessid"]) 
            if filename.endswith(".aln"):
                self.extract_info_from_fasta(filename=filename)
            find_seq3 = pd.merge(find_seq2, self.df, 
                         how = 'inner',
                         left_on = "accessid", 
                         right_on = "accessid")       
            find_seq3['aligned_sequence'] = find_seq3['aligned_sequence'].str.replace("-","")
            find_seq3 = find_seq3.rename(columns={"aligned_sequence": "sequence"})
            if i == 0:
                self.df_outlier_intra = find_seq3.copy()
            if i == 1:
                self.df_lower_outlier_inter = find_seq3.copy()
            if i == 2:
                self.df_upper_outlier_inter = find_seq3.copy()
        print("removing gaps from the sequences and saving an excel outlier file")
        self.df_outliers = [self.df_outlier_intra, self.df_lower_outlier_inter, self.df_upper_outlier_inter]
        df_outlier0 = pd.merge(self.df_outliers[0], self.df_outliers[1], 
                         how = 'left',
                         left_on = ["accessid","species","sequence"], 
                         right_on = ["accessid","species","sequence"])
        df_outlier = pd.merge(df_outlier0, self.df_outliers[2], 
                         how = 'left',
                         left_on = ["accessid","species","sequence"], 
                         right_on = ["accessid","species","sequence"])
        df_outlier = df_outlier.fillna(0)
        df_outlier["full_occurrence"] = df_outlier["occurrence_x"] + df_outlier["occurrence_y"] + df_outlier["occurrence"]
        df_outlier = df_outlier.sort_values('full_occurrence', ascending=False)
        self.df_outlier = df_outlier.reset_index(drop=True)
        if len(self.df_outlier)>1:
            print("tutti gli outliers si trovano in df_outlier")
            worst_outliers = df_outlier[df_outlier["full_occurrence"] > np.percentile(df_outlier["full_occurrence"],75)]
            self.worst_outliers = worst_outliers.reset_index(drop=True)
            if len(self.worst_outliers) < 5:
                self.worst_outliers = self.df_outlier
            if search_type == "full":
                for i in range(len(self.worst_outliers)):
                    stringa = ""
                    stringa = stringa + ">" 
                    stringa = stringa + self.worst_outliers["accessid"][i] + "_" + self.worst_outliers["species"][i]    
                    stringa = stringa + "\n" + self.worst_outliers["sequence"][i] 
                    file_fasta = stringa
                    dir = os.path.join(path,"TEMP")
                    if not os.path.exists(dir):
                        os.mkdir(dir)
                    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                    file_output = path + "TEMP//" + timestamp + "_" + str(i) + "_outliers.fasta"
                    fasta_output = file_output
                    o = open(file_output, "w")
                    o.write(stringa)
                    o.close
                    print(timestamp, " the fasta file ",i," is ready")
                if platform.system() == "Darwin" or platform.system() == "Linux":
                    onlyfiles = [f for f in listdir(path + "TEMP//") if isfile(join(path + "TEMP//", f))]
                if platform.system() == "Windows":
                    onlyfiles = [f for f in listdir(path + "TEMP\\") if isfile(join(path + "TEMP\\", f))]
                lista_fasta = [x for x in onlyfiles if x.endswith('outliers.fasta')] 
                self.lista_fasta = lista_fasta
                if platform.system() == "Darwin" or platform.system() == "Linux":
                    append_str = 'TEMP//'
                    lista_fasta = [append_str + sub for sub in lista_fasta]
                if platform.system() == "Windows":
                    append_str = 'TEMP\\'
                    lista_fasta = [append_str + sub for sub in lista_fasta]
                for i in range(len(lista_fasta)):   
                    f = open(path + lista_fasta[i], "r")
                    Query = f.read()
                    f.close
                    Query = Query.split("\n")[0]
                    print("Blasting ",Query)
                    self.execute_blast(query=lista_fasta[i], threshold=95, csv="none")
                    FromBlastDf = self.FromBlastDf[self.FromBlastDf[2] > 99]
                    if len(FromBlastDf) == 1:
                        FromBlastDf = self.FromBlastDf[self.FromBlastDf[2] > 95]
                        print("there are no species with ident >99")
                        os.remove(path+lista_fasta[i])
                        self.find_from_accessid(threshold="none", csv="none", gene="none", rank_type="s", filename=FromBlastDf[1], id_type="list", query="none")
                        print("The sequence ", Query, "has ident >95 and <99 with", list(set(self.df["species"])))
                    else:
                        os.remove(path+lista_fasta[i])
                        self.find_from_accessid(threshold="none", csv="none", gene="none", rank_type="s", filename=FromBlastDf[1], id_type="list", query="none")
                        print("The sequence ", Query, "has ident >99 with", list(set(self.df["species"])))
        else:
            print("no outlier found")
            
    def get_boxplot(self, title, rank):
        plt.boxplot((self.intra , self.inter)) 
        if rank == "Population":
            plt.xticks([1, 2], ['Intrapopulation', 'Interpopulation'])
        if rank == "Species":
            plt.xticks([1, 2], ['Intraspecific', 'Interspecific'])
        plt.title(title, fontsize=15)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        plt.savefig(path + "boxplot_" + str(title).replace(" ", "_") + timestamp + ".pdf", bbox_inches = 'tight')
        
    def barcoding_hist(self, comparison, rank, title):  
        comparison = comparison.split(",")
        loc = "upper left"
        y, x, _ = plt.hist(self.final_intra_range, bins = 1000) 
        intra_max = y.max()
        y, x, _ = plt.hist(self.final_inter_range, bins = 1000) 
        inter_max = y.max()
        intra_max, inter_max
        fig, ax = plt.subplots()
        if rank == "Species":
            ax.hist(self.final_inter_range, 1000, None, ec='grey', fc='grey', lw=1, 
                    histtype='step', label='Interspecific distances', fill = True)
            ax.hist(self.final_intra_range, 1000, None, ec='lightgrey', fc='lightgrey', lw=1, 
                    histtype='step', label='Intraspecific distances', fill = True)
        if rank == "Population":
            ax.hist(self.final_inter_range, 1000, None, ec='grey', fc='grey', lw=1, 
                    histtype='step', label='Interpopulation distances', fill = True)
            ax.hist(self.final_intra_range, 1000, None, ec='lightgrey', fc='lightgrey', lw=1, 
                    histtype='step', label='Intrapopulation distances', fill = True)
        #https://matplotlib.org/stable/gallery/color/named_colors.html
        ec = ['red','blue','orange','green','purple','aqua','indigo','yellow','violet','lightblue','hotpink','olive']
        fc = ec
        i = 0
        u=0
        m = 0
        for k in self.dict.keys():
            if rank == "Species":
                if (k.split("_")[0] in self.interest_specific_name) and (k.split("_")[1] 
                                                                         in self.interest_specific_name):
                    for c in range(len(comparison)): 

                        if(((comparison[c].split(" vs ")[0].split(" ")[1] == k.split("_")[0]) 
                            and (comparison[c].split(" vs ")[1].split(" ")[1] == k.split("_")[1])) 
                           or ((comparison[c].split(" vs ")[1].split(" ")[1] == k.split("_")[0]) 
                               and (comparison[c].split(" vs ")[0].split(" ")[1] == k.split("_")[1]))):
                            if len(self.dict[k])/inter_max > 0.04 and len(self.dict[k])/intra_max > 0.04: 
                                ax.hist(self.dict[k], 1000, None, ec=ec[i], fc=fc[i],  
                                        lw=1, histtype='step', label=comparison[c], fill = True)
                                i = i+1
                            else:
                                if max(self.dict[k]) - min(self.dict[k]) > max(self.final_intra_range):
                                    if min(self.dict[k]) == 0:
                                        plt.plot([min(self.dict[k])],[inter_max/80+u], marker = ("*"), c=ec[i], 
                                             label=comparison[c], linestyle="None")
                                        u= u+(inter_max/100)
                                    else:
                                        plt.plot([min(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                                 label=comparison[c], linestyle="None")
                                        plt.plot([max(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                                 label=comparison[c], linestyle="None")
                                else:
                                    plt.plot([mean(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                         label=comparison[c], linestyle="None")
                            i = i+1
            if rank == "Population":
                if (k.split("_")[0] in self.interest_specific_name) and (k.split("_")[1] 
                                                                     in self.interest_specific_name):
                    for c in range(len(comparison)): 
                        if(((comparison[c].split(" vs ")[0] == k.split("_")[0]) 
                            and (comparison[c].split(" vs ")[1] == k.split("_")[1])) 
                           or ((comparison[c].split(" vs ")[1] == k.split("_")[0]) 
                               and (comparison[c].split(" vs ")[0] == k.split("_")[1]))):
                            if len(self.dict[k])/inter_max > 0.04 and len(self.dict[k])/intra_max > 0.04: 
                                ax.hist(self.dict[k], 1000, None, ec=ec[i], fc=fc[i],  
                                        lw=1, histtype='step', label=comparison[c], fill = True)
                                i = i+1
                            else:
                                if max(self.dict[k]) - min(self.dict[k]) > max(self.final_intra_range):
                                    if min(self.dict[k]) == 0:
                                        plt.plot([min(self.dict[k])],[inter_max/80+u], marker = ("*"), c=ec[i], 
                                             label=comparison[c], linestyle="None")
                                        u= u+(inter_max/100)
                                    else:
                                        plt.plot([min(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                                 label=comparison[c], linestyle="None")
                                        plt.plot([max(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                                 label=comparison[c], linestyle="None")
                                else:
                                    plt.plot([mean(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                         label=comparison[c], linestyle="None")
                            i = i+1
                    
        '''ec = ['red','blue','orange','green','purple','aqua','indigo','yellow','violet','lightblue','hotpink','olive']
        fc = ec
        i = 0
        for k in self.dict.keys():
            if (k.split("_")[0] in self.interest_specific_name) and (k.split("_")[1] 
                                                                     in self.interest_specific_name):
                for c in range(len(comparison)): 
                    if(((comparison[c].split(" vs ")[0].split(" ")[1] == k.split("_")[0]) 
                        and (comparison[c].split(" vs ")[1].split(" ")[1] == k.split("_")[1])) 
                       or ((comparison[c].split(" vs ")[1].split(" ")[1] == k.split("_")[0]) 
                           and (comparison[c].split(" vs ")[0].split(" ")[1] == k.split("_")[1]))):
                        if len(self.dict[k])/inter_max > 0.04 and len(self.dict[k])/intra_max > 0.04: 
                            ax.hist(self.dict[k], 1000, None, ec=ec[i], fc=fc[i],  
                                    lw=1, histtype='step', label=comparison[c], fill = True)
                            i = i+1
                        else:
                            if max(self.dict[k]) - min(self.dict[k]) > max(self.final_intra_range):
                                plt.plot([min(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                         label=comparison[c], linestyle="None")
                                plt.plot([max(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                         label=comparison[c], linestyle="None")
                            else:
                                plt.plot([mean(self.dict[k])],[inter_max/80], marker = ("*"), c=ec[i], 
                                         label=comparison[c], linestyle="None")
                            i = i+1'''
        handles, labels = plt.gca().get_legend_handles_labels()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc=loc, prop={'size': 7}, frameon=False)
        ax.set_xlabel('Distance')
        ax.set_ylabel('Frequency')
        ax.set_title(title)
        plt.show()
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        fig.savefig(path + "Barcoding_Gap_" + str(title).replace(" ", "_") + timestamp + ".pdf", bbox_inches = 'tight')
    
    def Barcoding_analysis(self, filename, rank, kingdom, interesting_species, search_type, comparison, title):
        self.create_dataframe_matrix(filename)
        print("processing the pairwise distance matrix")
        self.process_matrix(rank, kingdom, interesting_species=interesting_species)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, "searching for the outliers")
        self.find_outlier(filename, search_type)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, "getting the boxplot")
        self.get_boxplot(title)
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        print(timestamp, "getting the barcoding histogram")
        self.barcoding_hist(comparison, rank, title)
 
    def iq_tree(self, filename):
        print("creating a Maximum Likehood tree..")
        subprocess.call(self.iqtree + " -s " + path + "//" + filename 
                        + " -st DNA -m TEST -bb 1000 -alrt 1000",shell=True)

