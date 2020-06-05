#!/usr/bin/python3
########################################################################
########################## traineeship, part 1: data collection using NCBI
##### Loading required package
### Environmental modules
import sys, os
import os.path
###
import re
# dnf install python3-mysql mysql-connector-python3
# pip3 install --user mysqlclient
import MySQLdb as my
# pip3 install PyMySQL
from eutils import Client
# pip3 install biopython
from Bio import Entrez
import csv
# pip3 install --user beautifulsoup4
import urllib.request
from bs4 import BeautifulSoup
######
print("\nPython version {}\n".format(sys.version))
print("\nThis script performs an online search in NCBI databases: RefSeq Nucleotide, dbVar, ClinVar & dbSNP.\nThe NCBI E-utilities (ESearch, ESummary) are used.\n")
print("RNA transcript variants for a specific gene will be gathered from UCSC (hg19-hg38) and NCBI RefSeq.\n")
gene = input("Specify the gene that will be used in the search query: ")
print("Choose one of the following options:\n")
print("1 ---> If you want to search for <Pathogenic Copy Number Variations in Humans (dbVar)>\n")
print("2 ---> If you want to search for <Copy Number Variations without clinical assertion in Humans (dbVar)>\n")
print("3 ---> If you want to search for <Copy Number Variations with clinical interpretation in Humans (dbVar)>\n")
print("4 ---> If you want to search for <Insertions in Humans (dbVar)>\n")
print("5 ---> If you want to search for <Inversions in Humans (dbVar)>\n")
print("6 ---> If you want to search for <Short Tandem Repeats in Humans (dbVar)>\n")
print("7 ---> If you want to search for <Less common types of variants in Humans (dbVar)>\n")
print("8 ---> If you want to search for <Short ClinVar variants>\n")
print("9 ---> If you want to search for <Long ClinVar variants>\n")
print("10 ---> If you want to search for <Common (>1%) short genetic variations in dbSNP>\n")
print("11 ---> If you want to search for <PubMed Cited short genetic variations in dbSNP>\n")
print("12 ---> If you want to search for 4, 5, 6, 7, 8, 9, 10, 11\n")
choiceofsearch = str(input("Number: "))
listofnumbers = ["1","2","3","4","5","6","7","8","9","10","11","12"]
### Check if number in listofnumbers
while choiceofsearch not in listofnumbers:
    choiceofsearch = str(input("Number mentioned above: "))
print("\nThe current Working Directory is '{}'.\nThe search results will be saved in the folder 'NCBI-search-results'.\n".format(os.getcwd()))
### Create new folder if it doesn't exist yet
newfolder = os.getcwd() + "/NCBI-search-results" 
newfolder2 = os.getcwd() + "/NCBI-search-results/" + gene
if os.path.isdir(newfolder):
    os.chdir(newfolder)
    if os.path.isdir(newfolder2):
        os.chdir(newfolder2)
    else:
        os.mkdir(newfolder2)
        os.chdir(newfolder2)
else:
    os.mkdir(newfolder)
    os.chdir(newfolder)
    os.mkdir(newfolder2)
    os.chdir(newfolder2)

#########################################################################################################

searchUCSC = {"hg19":'results-transcripts-UCSC-hg19.csv',"hg38":'results-transcripts-UCSC-hg38.csv'}
for key in searchUCSC:
    ##### Connection UCSC
    db = my.connect(host="genome-euro-mysql.soe.ucsc.edu",
        user="genomep",
        passwd="password",
        db=key)
    c = db.cursor()
    ##### ncbiRefSeq search
    no_rows = c.execute("SELECT * FROM ncbiRefSeq WHERE name2 LIKE '"+gene+"%'")
    result = c.fetchall()
    ##### Close database
    db.close()
    print("\nLoading currently available accession numbers from NCBI RefSeq table...")
    print("="*70)
    print("\nTranscript variant accession numbers: ")
    accList = []
    ### Save data to csv file
    with open(searchUCSC[key], mode='w') as result_transcripts:
        result_writer = csv.writer(result_transcripts,delimiter=';')
        result_writer.writerow(["chromosome","start","end","strand","gene","exon","transcript","symbol","ranges"])
        transcriptvar = 0
        for row in result:
            transcript = row[1]
            print(transcript)
            accList.append(transcript)
            starts = str(row[9])[2:-2]
            ends = str(row[10])[2:-2]
            starts1 = starts.split(",")
            ends1 = ends.split(",")
            j = 0
            ex = 1
            transcriptvar += 1
            for i in starts1:
                result_writer.writerow([row[2],i,ends1[j],row[3],row[12],ex,transcriptvar,row[1],str(row[4])+"-"+str(row[5])])
                j += 1
                ex += 1
    print("\nSearch results: {}\n".format(no_rows))
    ### Close csv file
    result_transcripts.close()
    print("Results are in '"+searchUCSC[key]+"'\n")

#########################################################################################################

##### API-key (NCBI)
eclient = Client(api_key="8ecce891e7fa036ff84bccc7c74e5138dc09")
#gene_efetch = eclient.efetch(db='gene', id=91039)
Entrez.email = "iris.raes@hotmail.com"

#########################################################################################################

##### nucleotide search
### Setting up query 
mRNAtranscripts = []
transcriptmRNA_esearch = eclient.esearch(db='nucleotide',
            term='('+gene+'[gene] AND "Homo sapiens"[Primary Organism] AND refseq[filter]) NOT biomol_genomic[PROP]')
print("\nLoading currently available ids from Entrez nucleotide...")
print("="*70)
print("\nTranscript variant ids: ")
print(transcriptmRNA_esearch.ids)
for item in transcriptmRNA_esearch.ids:
    mRNAtranscripts.append(item)
print("\nSearch results: {}\n".format(transcriptmRNA_esearch.count))
### Esummary for retrieving information
### For each id in mRNAtranscripts
### Save data to csv file
with open('results-nucleotide.csv', mode='w') as result_nucleotide:
    result_writer = csv.writer(result_nucleotide,delimiter=';')
    result_writer.writerow(["transcript id","description","transcript variant","accession","Chr","length in bp"])
    for ids in mRNAtranscripts:
        handle = Entrez.esummary(db="nucleotide", id=ids)
        record = Entrez.read(handle)
        handle.close()
        ### Write info to csv file, row by row
        splittedtitle = record[0]["Title"].split(",")
        result_writer.writerow([record[0]["Id"],splittedtitle[0],splittedtitle[1],record[0]["AccessionVersion"],"Chr19",str(record[0]["Length"])])
        ###
### Close csv file
result_nucleotide.close()
print("Results are in 'results-nucleotide.csv'\n")

#########################################################################################################

if choiceofsearch == "12":
    terms = {gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "insertion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-insertion-dbVar.csv',
    gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "inversion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-inversion-dbVar.csv',
    gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "short tandem repeat"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-STR-dbVar.csv',
    gene+'[All Fields] AND ("Homo sapiens"[Organism] NOT "copy number variation"[Variant Type] NOT "insertion"[Variant Type] NOT "short tandem repeat"[Variant Type] NOT "inversion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-lesscommon-dbVar.csv'}
if choiceofsearch == "1":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] AND "Pathogenic"[clinical_interpretation]) AND "VARIANT"[OBJ_TYPE]':'results-CNV-dbVar.csv'}
if choiceofsearch == "2":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] AND "not reported"[clinical_interpretation]) AND "VARIANT"[OBJ_TYPE]':'results-CNV-noclinint-dbVar.csv'}
if choiceofsearch == "3":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] NOT "not reported"[clinical_interpretation]) AND "VARIANT"[OBJ_TYPE]':'results-CNV-clinint-dbVar.csv'}
if choiceofsearch == "4":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "insertion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-insertion-dbVar.csv'}
if choiceofsearch == "5":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "inversion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-inversion-dbVar.csv'}
if choiceofsearch == "6":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "short tandem repeat"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-STR-dbVar.csv'}
if choiceofsearch == "7":
    terms={gene+'[All Fields] AND ("Homo sapiens"[Organism] NOT "copy number variation"[Variant Type] NOT "insertion"[Variant Type] NOT "short tandem repeat"[Variant Type] NOT "inversion"[Variant Type]) AND "VARIANT"[OBJ_TYPE]':'results-lesscommon-dbVar.csv'}

##### dbVar search
### Setting up query
if choiceofsearch in ["12","1","2","3","4","5","6","7"]:
    for key in terms:
        dbVar = []
        dbVar_esearch = eclient.esearch(db='dbVar',term=key)
        print("\nLoading currently available ids from dbVar...")
        print("="*70)
        print("dbVar ids: ")
        print(dbVar_esearch.ids)
        for item in dbVar_esearch.ids:
            dbVar.append(item)
        print("\nSearch results: {}\n".format(dbVar_esearch.count))
        ### Esummary for retrieving information
        ### For each id in dbVar
        ### Save data to csv file
        with open(terms[key], mode='w') as result_dbVar:
            result_writer = csv.writer(result_dbVar,delimiter=';')
            if choiceofsearch in ["1","3"]:
                result_writer.writerow(["dbVar variant id","variant region id","type","variant call id","variant call type","copy number","allele origin","subject phenotype","study ID","clinical assertion","ClinVar id","Chr","assembly1","assembly2","assembly3"])
            else:
                result_writer.writerow(["dbVar variant id","variant region id","type","variant call id","variant call type","sample id","zygosity","subject phenotype","study ID","clinical assertion","Chr","assembly1","assembly2","assembly3"])
            for ids in dbVar:
                handle = Entrez.esummary(db="dbVar", id=ids)
                record = Entrez.read(handle)
                handle.close()
                varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
                if choiceofsearch in ["1","3"]:
                    url = "https://www.ncbi.nlm.nih.gov/dbvar/variants/"+varregid+"/#VariantClinical"
                    response = urllib.request.urlopen(url)
                    content = response.read()
                    soup = BeautifulSoup(content, 'html.parser')
                    table = soup.find('table', {"data-column-types" : '["fnc:strnumSort","fnc:strnumSort"]'})
                    headers = table.find('thead')
                    values = table.find('tbody')
                    header = headers.find_all('th')
                    value = values.find_all('td')
                    value = [re.sub('<a.*?>|</a>', '', str(v)) for v in value]
                    value = [re.sub('<td.*?>|</td>', '', str(v)) for v in value]
                    if len(header) < len(value):
                        val = value[0:len(header)*2]
                        value1 = val[:len(val)//2]
                        value2 = val[len(val)//2:]
                        if len(value)/len(header) > 2:
                            merged = [v1+"//"+v2+"..." for v1,v2 in zip(value1,value2)]
                        else:
                            merged = [v1+"//"+v2 for v1,v2 in zip(value1,value2)]
                        merged = [re.sub(',', '-', str(v)) for v in merged]
                        dictionary = {h.string:v for h,v in zip(header,merged)}
                    else:
                        dictionary = {h.string:v for h,v in zip(header,value)}
                    varcall = dictionary.get("Variant Call ID","")
                    varcalltype = dictionary.get("Type","")
                    allele = dictionary.get("Allele Origin","")
                    phen = dictionary.get("Subject Phenotype","")
                    clinvarid = dictionary.get("ClinVar ID","")
                    copynr = dictionary.get("Copy number","")
                else:
                    url = "https://www.ncbi.nlm.nih.gov/dbvar/variants/"+varregid+"/#VariantDetails"
                    response = urllib.request.urlopen(url)
                    content = response.read()
                    soup = BeautifulSoup(content, 'html.parser')
                    table = soup.find('div', {"id" : "support_variant_set"})
                    headers = table.find('thead')
                    values = table.find('tbody')
                    header = headers.find_all('th')
                    value = values.find_all('td')
                    value = [re.sub('<a.*?>|</a>', '', str(v)) for v in value]
                    value = [re.sub('<td.*?>|</td>', '', str(v)) for v in value]
                    if len(header) < len(value):
                        val = value[0:len(header)*2]
                        value1 = val[:len(val)//2]
                        value2 = val[len(val)//2:]
                        if len(value)/len(header) > 2:
                            merged = [v1+"//"+v2+"..." for v1,v2 in zip(value1,value2)]
                        else:
                            merged = [v1+"//"+v2 for v1,v2 in zip(value1,value2)]
                        merged = [re.sub(',', '-', str(v)) for v in merged]
                        dictionary = {h.string:v for h,v in zip(header,merged)}
                    else:
                        dictionary = {h.string:v for h,v in zip(header,value)}
                    varcall = dictionary.get("Variant Call ID","")
                    varcalltype = dictionary.get("Type","")
                    sample = dictionary.get("Sample ID","")
                    zygosity = dictionary.get("Zygosity","")
                    phen = dictionary.get("Subject Phenotype","")
                types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
                studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
                try:
                    clinicalassertion = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarClinicalSignificanceList')
                except:
                    clinicalassertion = "n/a"
                Chr = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
                try:
                    assembly3 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][2].get('Assembly')
                    start3 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][2].get('Chr_start')
                    end3 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][2].get('Chr_end')
                except:
                    assembly3 = "not applicable"
                    start3 = "X"
                    end3 = "X"
                ### Write info to csv file, row by row
                if choiceofsearch in ["1","3"]:
                    result_writer.writerow([ids,varregid,types,varcall,varcalltype,copynr,allele,phen,studyid,
                    clinicalassertion,clinvarid,Chr,assembly1+":"+start1+"-"+end1,assembly2+":"+start2+"-"+end2,assembly3+":"+start3+"-"+end3])
                else:
                    result_writer.writerow([ids,varregid,types,varcall,varcalltype,sample,zygosity,phen,studyid,
                    clinicalassertion,Chr,assembly1+":"+start1+"-"+end1,assembly2+":"+start2+"-"+end2,assembly3+":"+start3+"-"+end3])
        ### Close csv file
        result_dbVar.close()
        print("Results are in '"+terms[key]+"'\n")

#########################################################################################################

if choiceofsearch == "12":
    terms = {gene+'[gene] AND "Single nucleotide"':'results-ClinVar-short.csv',gene+'[gene] NOT "Single nucleotide"':'results-ClinVar-long.csv'}
if choiceofsearch == "8":
    terms={gene+'[gene] AND "Single nucleotide"':'results-ClinVar-short.csv'}
if choiceofsearch == "9":
    terms={gene+'[gene] NOT "Single nucleotide"':'results-ClinVar-long.csv'}

##### ClinVar search
### Setting up query 
if choiceofsearch in ["12","8","9"]:
    for key in terms:
        ClinVar = []
        ClinVar_esearch = eclient.esearch(db='ClinVar',term=key)
        print("\nLoading currently available ids from ClinVar...")
        print("="*70)
        print("\nClinVar ids: ")
        print(ClinVar_esearch.ids)
        for item in ClinVar_esearch.ids:
            ClinVar.append(item)
        print("\nSearch results: {}\n".format(ClinVar_esearch.count))
        ### Esummary for retrieving information
        ### For each id in ClinVar
        ### Save data to csv file
        with open(terms[key], mode='w') as result_ClinVar:
            result_writer = csv.writer(result_ClinVar,delimiter=';')
            result_writer.writerow(["ClinVar variant id","title","accession","description","type","protein change","Chr","assembly1","assembly2","assembly3","source id"])
            for ids in ClinVar:
                handle = Entrez.esummary(db="ClinVar", id=ids)
                record = Entrez.read(handle)
                handle.close()
                title = record['DocumentSummarySet']['DocumentSummary'][0].get('title')
                accession = record['DocumentSummarySet']['DocumentSummary'][0].get('accession')
                types = record['DocumentSummarySet']['DocumentSummary'][0].get('obj_type')
                description = record['DocumentSummarySet']['DocumentSummary'][0]['clinical_significance'].get('description')
                protein_change = record['DocumentSummarySet']['DocumentSummary'][0].get('protein_change')
                if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'] != []:
                    Chr = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('chr')
                    assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('assembly_name')
                    start1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('start')
                    end1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('stop')
                    if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0].get('variation_loc') != []:
                        try:
                            assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('assembly_name')
                            start2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('start')
                            end2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('stop')
                        except:
                            assembly2 = "not applicable"
                            start2 = "X"
                            end2 = "X"
                        try:
                            assembly3 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][2].get('assembly_name')
                            start3 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][2].get('start')
                            end3 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][2].get('stop')
                        except:
                            assembly3 = "not applicable"
                            start3 = "X"
                            end3 = "X"
                if record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'] != []:
                    try:
                        dbsource = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0].get('db_source')
                        dbid = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0].get('db_id')
                        dbid2 = ""
                    except:
                        dbsource = ""
                        dbid = ""
                if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_xrefs'] != []:
                    try:
                        dbsource = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_xrefs'][0].get('db_source')
                        dbid = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_xrefs'][0].get('db_id')
                        dbid2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_xrefs'][1].get('db_id')
                    except:
                        dbsource = ""
                        dbid = ""
                        dbid2 = ""
                ### Write info to csv file, row by row
                result_writer.writerow([ids,title,accession,description,types,protein_change,Chr,assembly1+":"+start1+"-"+end1,assembly2+":"+start2+"-"+end2,assembly3+":"+start3+"-"+end3,dbsource+" ("+dbid+" "+dbid2+")"])
                ###
        ### Close csv file
        result_ClinVar.close()
        print("Results are in '"+terms[key]+"'\n")
    
#########################################################################################################

if choiceofsearch  == "12":
    terms = {gene+'[All Fields] AND (00000.0100[GLOBAL_MAF] : 00000.1000[GLOBAL_MAF])':'results-common-dbSNP.csv',gene+'[All Fields] AND snp_pubmed_cited[Filter]':'results-cited-dbSNP.csv'}
if choiceofsearch == "10":
    terms = {gene+'[All Fields] AND (00000.0100[GLOBAL_MAF] : 00000.1000[GLOBAL_MAF])':'results-common-dbSNP.csv'}
if choiceofsearch == "11":
    terms = {gene+'[All Fields] AND snp_pubmed_cited[Filter]':'results-cited-dbSNP.csv'}

##### dbSNP search
### Setting up query 
if choiceofsearch in ["12","10","11"]:
    for key in terms:
        dbSNP = []
        dbSNP_esearch = eclient.esearch(db='snp',term=key)
        print("\nLoading currently available ids from dbSNP...")
        print("="*70)
        print("\ndbSNP ids: ")
        print(dbSNP_esearch.ids)
        for item in dbSNP_esearch.ids:
            dbSNP.append(item)
        print("\nSearch results: {}\n".format(dbSNP_esearch.count))
        ### Esummary for retrieving information
        ### For each id in dbSNP
        ### Save data to csv file
        with open(terms[key], mode='w') as result_dbSNP:
            result_writer = csv.writer(result_dbSNP,delimiter=';')
            result_writer.writerow(["dbSNP variant id","allele","variation type","functional consequence","Chr","assembly1","assembly2","PMID","title","author"])
            for ids in dbSNP:
                handle = Entrez.esummary(db="snp", id=ids)
                record = Entrez.read(handle)
                handle.close()
                docsum = record['DocumentSummarySet']['DocumentSummary'][0].get('DOCSUM')
                allele =  re.search(r"SEQ=\[.*?\]", docsum).group(0)
                vartype = record['DocumentSummarySet']['DocumentSummary'][0].get('SNP_CLASS')
                conseq = record['DocumentSummarySet']['DocumentSummary'][0].get('FXN_CLASS')
                conseq = conseq.replace(",", "/")
                Chr = record['DocumentSummarySet']['DocumentSummary'][0].get('CHR')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0].get('CHRPOS_PREV_ASSM')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0].get('CHRPOS')
                try:
                    url = "https://www.ncbi.nlm.nih.gov/snp/"+str(ids)+"#publications"
                    response = urllib.request.urlopen(url)
                    content = response.read()
                    soup = BeautifulSoup(content, 'html.parser')
                    table = soup.find('table', {"id" : 'publication_datatable'})
                    headers = table.find('thead')
                    values = table.find('tbody')
                    header = headers.find_all('th')
                    value = values.find_all('td')
                    value = [re.sub('<a.*?>|</a>', '', str(v)) for v in value]
                    value = [re.sub('<td.*?>|</td>', '', str(v)) for v in value]
                    value = [re.sub('\n', '', str(v)) for v in value]
                    if len(header) < len(value):
                        val = value[0:len(header)*2]
                        value1 = val[:len(val)//2]
                        value2 = val[len(val)//2:]
                        if len(value)/len(header) > 2:
                            merged = [v1+"//"+v2+"..." for v1,v2 in zip(value1,value2)]
                        else:
                            merged = [v1+"//"+v2 for v1,v2 in zip(value1,value2)]
                        merged = [re.sub(',', '-', str(v)) for v in merged]
                        dictionary = {h.string:v for h,v in zip(header,merged)}
                    else:
                        dictionary = {h.string:v for h,v in zip(header,value)}
                    pmid = dictionary.get("PMID","")
                    title = dictionary.get("Title","")
                    author = dictionary.get("Author","")
                except:
                    pmid = ""
                    title = ""
                    author = ""
                ### Write info to csv file, row by row
                result_writer.writerow(["rs"+str(ids),allele,vartype,conseq,Chr,"GRCh37-"+assembly1,"GRCh38-"+assembly2,pmid,title,author])
                ###
        ### Close csv file
        result_dbSNP.close()
        print("Results are in '"+terms[key]+"'\n")

print("\n\n\t*** NCBI Search successful ***\n\n")
