#!/usr/bin/python3
########################################################################
########################## traineeship, part 1: data collection using NCBI
##### Loading required package
### Environmental modules
import sys, os
import os.path
###
# dnf install python3-mysql mysql-connector-python3
# pip3 install --user mysqlclient
import MySQLdb as my
# pip3 install PyMySQL
from eutils import Client
# pip3 install biopython
from Bio import Entrez
import csv
######
print("\nPython version {}\n".format(sys.version))
print("\nThis script performs an online search in NCBI databases: RefSeq, Nucleotide, dbVar & ClinVar.\nThe NCBI E-utilities (ESearch, ESummary) are used.\n")
print("RNA transcript variants for a specific gene will be gathered from UCSC (hg19) and NCBI RefSeq.\n")
gene = input("Specify the gene that will be used in the search query: ")
print("Choose one of the following options:\n")
print("1 ---> If you want to search for <Pathogenic Copy Number Variations in Human (dbVar)>\n")
print("2 ---> If you want to search for <Copy Number Variations without clinical assertion in Human (dbVar)>\n")
print("3 ---> If you want to search for <Copy Number Variations not reported as 'Pathogenic' in Human (dbVar)>\n")
print("4 ---> If you want to search for <Insertions in Human (dbVar)>\n")
print("5 ---> If you want to search for <Inversions in Human (dbVar)>\n")
print("6 ---> If you want to search for <Less common insertions and deletions in Human (dbVar)>\n")
print("7 ---> If you want to search for <Short Tandem Repeats in Human (dbVar)>\n")
print("8 ---> If you want to search for <Substitutions, alterations, dublications and translocations in Human (dbVar)>\n")
print("9 ---> If you want to search for <Short ClinVar variants>\n")
print("10 ---> If you want to search for <Long ClinVar variants>\n")
print("11 ---> If you want to search for 1, 2, 3, 4, 5, 6, 7, 8, 9, 10\n")
choiceofsearch = str(input("Number: "))
listofnumbers = ["1","2","3","4","5","6","7","8","9","10","11"]
while choiceofsearch not in listofnumbers:
    choiceofsearch = str(input("Number mentioned above: "))
print("\nThe current Working Directory is '{}'.\nThe search results will be saved in the folder 'NCBI-search-results'.\n".format(os.getcwd()))
newfolder = os.getcwd() + "/NCBI-search-results"
if os.path.isdir(newfolder):
    os.chdir(newfolder)
else:
    os.mkdir(newfolder)
    os.chdir(newfolder)
######
##### Connection UCSC
db = my.connect(host="genome-euro-mysql.soe.ucsc.edu",
   user="genomep",
   passwd="password",
   db="hg19")
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
with open('results-transcripts-UCSC.csv', mode='w') as result_transcripts:
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
print("Results are in 'results-transcripts-UCSC.csv'\n")

##### API-key (NCBI)
eclient = Client(api_key="8ecce891e7fa036ff84bccc7c74e5138dc09")
#gene_efetch = eclient.efetch(db='gene', id=91039)

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
Entrez.email = "iris.raes@hotmail.com"
### For each id in mRNAtranscripts
### Save data to csv file
with open('results-nucleotide.csv', mode='w') as result_nucleotide:
    result_writer = csv.writer(result_nucleotide,delimiter=';')
    result_writer.writerow(["transcript_id","description","transcript_variant","accession","Chr","length_in_bp"])
    for ids in mRNAtranscripts:
        handle = Entrez.esummary(db="nucleotide", id=ids)
        record = Entrez.read(handle)
        handle.close()
        ### Write info to csv file, row by row
        splittedtitle = record[0]["Title"].split(",")
        result_writer.writerow([record[0]["Id"],splittedtitle[0],splittedtitle[1],record[0]["AccessionVersion"],"Chr19","1-"+str(record[0]["Length"])])
        ###
### Close csv file
result_nucleotide.close()
print("Results are in 'results-nucleotide.csv'\n")

if choiceofsearch == "1" or choiceofsearch == "11":
    ##### dbVar search
    ### Setting up query
    ### Pathogenic Copy Number Variation in human
    CNV = []
    CNV_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] AND "Pathogenic"[clinical_interpretation])')
    print("\nLoading currently available ids from dbVar...")
    print("="*70)
    print("dbVar ids: ")
    print(CNV_esearch.ids)
    for item in CNV_esearch.ids:
        CNV.append(item)
    print("\nSearch results: {}\n".format(CNV_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in CNV
    ### Save data to csv file
    with open('results-CNV-dbVar.csv', mode='w') as result_CNV:
        result_writer = csv.writer(result_CNV,delimiter=';')
        result_writer.writerow(["CNV_variant_id","variant_region_id","type","study_ID","clinical_assertion","Chr_1","assembly1","Chr_2","assembly2"])
        for ids in CNV:
            handle = Entrez.esummary(db="dbVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
            studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
            clinicalassertion = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarClinicalSignificanceList')
            if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
            ### Write info to csv file, row by row
            result_writer.writerow([ids,varregid,types,studyid,clinicalassertion,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
            ###
    ### Close csv file
    result_CNV.close()
    print("Results are in 'results-CNV-dbVar.csv'\n")

if choiceofsearch == "2" or choiceofsearch == "11":
    ##### dbVar search
    ### Setting up query
    ### Copy Number Variation without clinical assertion in human
    CNV = []
    CNV_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] AND "not reported"[clinical_interpretation])')
    print("\nLoading currently available ids from dbVar...")
    print("="*70)
    print("dbVar ids: ")
    print(CNV_esearch.ids)
    for item in CNV_esearch.ids:
        CNV.append(item)
    print("\nSearch results: {}\n".format(CNV_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in CNV
    ### Save data to csv file
    with open('results-CNV-noclin-dbVar.csv', mode='w') as result_CNV:
        result_writer = csv.writer(result_CNV,delimiter=';')
        result_writer.writerow(["CNV_variant_id","variant_region_id","type","study_ID","clinical_assertion","Chr_1","assembly1","Chr_2","assembly2"])
        for ids in CNV:
            handle = Entrez.esummary(db="dbVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
            studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
            clinicalassertion = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarClinicalSignificanceList')
            if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
            ### Write info to csv file, row by row
            result_writer.writerow([ids,varregid,types,studyid,clinicalassertion,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
            ###
    ### Close csv file
    result_CNV.close()
    print("Results are in 'results-CNV-noclin-dbVar.csv'\n")

if choiceofsearch == "3" or choiceofsearch == "11":
    ##### dbVar search
    ### Setting up query
    ### Copy Number Variation not reported as pathogenic in human
    CNV = []
    listofCNV = ["Likely pathogenic","Uncertain significance","Benign","Likely benign","Benign/Likely benign","Conflicting interpretations of pathogenicity",
                        "Pathogenic/Likely pathogenic","association","conflicting data from submitters","drug response","protective","risk factor"]
    for searchitem in listofCNV:
        CNV_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "copy number variation"[Variant Type] AND "'+searchitem+'"[clinical_interpretation])')
        print("\nLoading currently available ids for <"+searchitem+"> from dbVar...")
        print("="*70)
        print("dbVar ids: ")
        print(CNV_esearch.ids)
        for item in CNV_esearch.ids:
            CNV.append(item)
        print("\nSearch results: {}\n".format(CNV_esearch.count))
        ### Esummary for retrieving information
        Entrez.email = "iris.raes@hotmail.com"
        ### For each id in CNV
        ### Save data to csv file
        with open('results-CNV-notpathogenic-dbVar.csv', mode='w') as result_CNV:
            result_writer = csv.writer(result_CNV,delimiter=';')
            result_writer.writerow(["CNV_variant_id","variant_region_id","type","study_ID","clinical_assertion","Chr_1","assembly1","Chr_2","assembly2"])
            for ids in CNV:
                handle = Entrez.esummary(db="dbVar", id=ids)
                record = Entrez.read(handle)
                handle.close()
                varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
                types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
                studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
                clinicalassertion = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarClinicalSignificanceList')
                if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                    Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                    assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                    start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                    end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                    Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                    assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                    start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                    end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
                ### Write info to csv file, row by row
                result_writer.writerow([ids,varregid,types,studyid,clinicalassertion,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
                ###
        ### Close csv file
        result_CNV.close()
    print("Results are in 'results-CNV-notpathogenic-dbVar.csv'\n")

if choiceofsearch == "4" or choiceofsearch == "11":
    ### Insertion in human
    insertion = []
    insertion_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "insertion"[Variant Type])')
    print("\nLoading currently available ids from dbVar...")
    print("="*70)
    print("dbVar ids: ")
    print(insertion_esearch.ids)
    for item in insertion_esearch.ids:
        insertion.append(item)
    print("\nSearch results: {}\n".format(insertion_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in insertion
    ### Save data to csv file
    with open('results-insertion-dbVar.csv', mode='w') as result_insertion:
        result_writer = csv.writer(result_insertion,delimiter=';')
        result_writer.writerow(["insertion_variant_id","variant_region_id","type","study_ID","Chr_1","assembly1","Chr_2","assembly2"])
        for ids in insertion:
            handle = Entrez.esummary(db="dbVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
            studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
            if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
            ### Write info to csv file, row by row
            result_writer.writerow([ids,varregid,types,studyid,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
            ###
    ### Close csv file
    result_insertion.close()
    print("Results are in 'results-insertion-dbVar.csv'\n")

if choiceofsearch == "5" or choiceofsearch == "11":
    ### Inversion in human
    inversion = []
    inversion_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "inversion"[Variant Type])')
    print("\nLoading currently available ids from dbVar...")
    print("="*70)
    print("dbVar ids: ")
    print(inversion_esearch.ids)
    for item in inversion_esearch.ids:
        inversion.append(item)
    print("\nSearch results: {}\n".format(inversion_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in inversion
    ### Save data to csv file
    with open('results-inversion-dbVar.csv', mode='w') as result_inversion:
        result_writer = csv.writer(result_inversion,delimiter=';')
        result_writer.writerow(["inversion_variant_id","variant_region_id","type","study_ID","Chr_1","assembly1","Chr_2","assembly2"])
        for ids in inversion:
            handle = Entrez.esummary(db="dbVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
            studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
            if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarRemappedAssemblyList'] != []:
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarRemappedAssemblyList'][0]
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarRemappedAssemblyList'][0]
            ### Write info to csv file, row by row
            result_writer.writerow([ids,varregid,types,studyid,"Chr19",assembly1,"Chr19",assembly2])
            ###
    ### Close csv file
    result_inversion.close()
    print("Results are in 'results-inversion-dbVar.csv'\n")

if choiceofsearch == "6" or choiceofsearch == "11":
    ##### dbVar search
    ### Setting up query
    ### Less common insertions and deletions in human
    indel = []
    listofindel = ["mobile element insertion","novel sequence insertion","mobile element deletion","delins"]
    for searchitem in listofindel:
        indel_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "'+searchitem+'"[Variant Type])')
        print("\nLoading currently available ids for <"+searchitem+"> from dbVar...")
        print("="*70)
        print("dbVar ids: ")
        print(indel_esearch.ids)
        for item in indel_esearch.ids:
            indel.append(item)
        print("\nSearch results: {}\n".format(indel_esearch.count))
        ### Esummary for retrieving information
        Entrez.email = "iris.raes@hotmail.com"
        ### For each id in indel
        ### Save data to csv file
        with open('results-lesscommon-indel-dbVar.csv', mode='w') as result_indel:
            result_writer = csv.writer(result_indel,delimiter=';')
            result_writer.writerow(["indel_variant_id","variant_region_id","type","study_ID","Chr_1","assembly1","Chr_2","assembly2"])
            for ids in indel:
                handle = Entrez.esummary(db="dbVar", id=ids)
                record = Entrez.read(handle)
                handle.close()
                varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
                types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
                studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
                if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                    Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                    assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                    start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                    end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                    Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                    assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                    start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                    end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
                else:
                    Chr_1 = ""
                    assembly1 = "not applicable"
                    start1 = "X"
                    end1 = "X"
                    Chr_2 = ""
                    assembly2 = "not applicable"
                    start2 = "X"
                    end2 = "X"
                ### Write info to csv file, row by row
                result_writer.writerow([ids,varregid,types,studyid,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
                ###
        ### Close csv file
        result_indel.close()
    print("Results are in 'results-lesscommon-indel-dbVar.csv'\n")

if choiceofsearch == "7" or choiceofsearch == "11":
    ### Short Tandem Repeats in human (seems to be less important)
    STR = []
    STR_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "short tandem repeat"[Variant Type])')
    print("\nLoading currently available ids from dbVar...")
    print("="*70)
    print("dbVar ids: ")
    print(STR_esearch.ids)
    for item in STR_esearch.ids:
        STR.append(item)
    print("\nSearch results: {}\n".format(STR_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in STR
    ### Save data to csv file
    with open('results-STR-dbVar.csv', mode='w') as result_STR:
        result_writer = csv.writer(result_STR,delimiter=';')
        result_writer.writerow(["STR_variant_id","variant_region_id","type","study_ID","Chr_1","assembly1","Chr_2","assembly2"])
        for ids in STR:
            handle = Entrez.esummary(db="dbVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
            studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
            if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
            ### Write info to csv file, row by row
            result_writer.writerow([ids,varregid,types,studyid,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
            ###
    ### Close csv file
    result_STR.close()
    print("Results are in 'results-STR-dbVar.csv'\n")

if choiceofsearch == "8" or choiceofsearch == "11":
    ##### dbVar search
    ### Setting up query
    ### Substitutions, alterations, dublications and translocations in human
    complexalt = []
    listofcomplexalt = ["complex substitution","sequence alteration","tandem duplication","translocation","complex chromosomal rearrangement"]
    for searchitem in listofcomplexalt:
        complexalt_esearch = eclient.esearch(db='dbVar',
            term=gene+'[All Fields] AND ("Homo sapiens"[Organism] AND "'+searchitem+'"[Variant Type])')
        print("\nLoading currently available ids for <"+searchitem+"> from dbVar...")
        print("="*70)
        print("dbVar ids: ")
        print(complexalt_esearch.ids)
        for item in complexalt_esearch.ids:
            complexalt.append(item)
        print("\nSearch results: {}\n".format(complexalt_esearch.count))
        ### Esummary for retrieving information
        Entrez.email = "iris.raes@hotmail.com"
        ### For each id in complexalt
        ### Save data to csv file
        with open('results-complexalt-dbVar.csv', mode='w') as result_complexalt:
            result_writer = csv.writer(result_complexalt,delimiter=';')
            result_writer.writerow(["indel_variant_id","variant_region_id","type","study_ID","Chr_1","assembly1","Chr_2","assembly2"])
            for ids in complexalt:
                handle = Entrez.esummary(db="dbVar", id=ids)
                record = Entrez.read(handle)
                handle.close()
                varregid = record['DocumentSummarySet']['DocumentSummary'][0].get('SV')
                types = record['DocumentSummarySet']['DocumentSummary'][0].get('dbVarVariantTypeList')
                studyid = record['DocumentSummarySet']['DocumentSummary'][0].get('ST')
                if record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'] != []:
                    Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr')
                    assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Assembly')
                    start1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_start')
                    end1 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][0].get('Chr_end')
                    Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr')
                    assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Assembly')
                    start2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_start')
                    end2 = record['DocumentSummarySet']['DocumentSummary'][0]['dbVarPlacementList'][1].get('Chr_end')
                else:
                    Chr_1 = ""
                    assembly1 = "not applicable"
                    start1 = "X"
                    end1 = "X"
                    Chr_2 = ""
                    assembly2 = "not applicable"
                    start2 = "X"
                    end2 = "X"
                ### Write info to csv file, row by row
                result_writer.writerow([ids,varregid,types,studyid,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2])
                ###
        ### Close csv file
        result_complexalt.close()
    print("Results are in 'results-complexalt-dbVar.csv'\n")

if choiceofsearch == "9" or choiceofsearch == "11":
    ##### ClinVar search
    ### Setting up query 
    ClinVar = []
    ClinVar_esearch = eclient.esearch(db='ClinVar',
            term=gene+'[gene] AND "Single gene"')
    print("\nLoading currently available ids from ClinVar...")
    print("="*70)
    print("\nClinVar ids: ")
    print(ClinVar_esearch.ids)
    for item in ClinVar_esearch.ids:
        ClinVar.append(item)
    print("\nSearch results: {}\n".format(ClinVar_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in ClinVar
    ### Save data to csv file
    with open('results-ClinVar-short.csv', mode='w') as result_ClinVar:
        result_writer = csv.writer(result_ClinVar,delimiter=';')
        result_writer.writerow(["ClinVar_variant_id","title","accession","type","description","protein_change","Chr_1","assembly1","Chr_2","assembly2","source_id"])
        for ids in ClinVar:
            handle = Entrez.esummary(db="ClinVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            title = record['DocumentSummarySet']['DocumentSummary'][0].get('title')
            accession = record['DocumentSummarySet']['DocumentSummary'][0].get('accession_version')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('obj_type')
            description = record['DocumentSummarySet']['DocumentSummary'][0]['clinical_significance'].get('description')
            protein_change = record['DocumentSummarySet']['DocumentSummary'][0].get('protein_change')
            if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('assembly_name')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('stop')
                if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0].get('variation_loc') != []:
                    try:
                        Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('chr')
                        assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('assembly_name')
                        start2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('start')
                        end2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('stop')
                    except:
                        Chr_2 = ""
                        assembly2 = "not applicable"
                        start2 = "X"
                        end2 = "X"
            if record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'] != []:
                dbsource = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0]['db_source']
                dbid = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0]['db_id']
            ### Write info to csv file, row by row
            result_writer.writerow([ids,title,accession,types,description,protein_change,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2,dbsource+" ("+dbid+")"])
            ###
    ### Close csv file
    result_ClinVar.close()
    print("Results are in 'results-ClinVar-short.csv'\n")

if choiceofsearch == "10" or choiceofsearch == "11":
    ##### ClinVar search
    ### Setting up query 
    ClinVar = []
    ClinVar_esearch = eclient.esearch(db='ClinVar',
            term=gene+'[gene] NOT "Single gene"')
    print("\nLoading currently available ids from ClinVar...")
    print("="*70)
    print("\nClinVar ids: ")
    print(ClinVar_esearch.ids)
    for item in ClinVar_esearch.ids:
        ClinVar.append(item)
    print("\nSearch results: {}\n".format(ClinVar_esearch.count))
    ### Esummary for retrieving information
    Entrez.email = "iris.raes@hotmail.com"
    ### For each id in ClinVar
    ### Save data to csv file
    with open('results-ClinVar-long.csv', mode='w') as result_ClinVar:
        result_writer = csv.writer(result_ClinVar,delimiter=';')
        result_writer.writerow(["ClinVar_variant_id","title","accession","type","description","protein_change","Chr_1","assembly1","Chr_2","assembly2","source_id"])
        for ids in ClinVar:
            handle = Entrez.esummary(db="ClinVar", id=ids)
            record = Entrez.read(handle)
            handle.close()
            title = record['DocumentSummarySet']['DocumentSummary'][0].get('title')
            accession = record['DocumentSummarySet']['DocumentSummary'][0].get('accession_version')
            types = record['DocumentSummarySet']['DocumentSummary'][0].get('obj_type')
            description = record['DocumentSummarySet']['DocumentSummary'][0]['clinical_significance'].get('description')
            protein_change = record['DocumentSummarySet']['DocumentSummary'][0].get('protein_change')
            if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'] != []:
                Chr_1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('chr')
                assembly1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('assembly_name')
                start1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('start')
                end1 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][0].get('stop')
                if record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0].get('variation_loc') != []:
                    try:
                        Chr_2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('chr')
                        assembly2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('assembly_name')
                        start2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('start')
                        end2 = record['DocumentSummarySet']['DocumentSummary'][0]['variation_set'][0]['variation_loc'][1].get('stop')
                    except:
                        Chr_2 = ""
                        assembly2 = "not applicable"
                        start2 = "X"
                        end2 = "X"
            if record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'] != []:
                try:
                    dbsource = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0].get('db_source')
                    dbid = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][0]['trait_xrefs'][0].get('db_id')
                except:
                    dbsource = "/"
                    dbid = ""
            ### Write info to csv file, row by row
            result_writer.writerow([ids,title,accession,types,description,protein_change,Chr_1,assembly1+":"+start1+"-"+end1,Chr_2,assembly2+":"+start2+"-"+end2,dbsource+" ("+dbid+")"])
            ###
    ### Close csv file
    result_ClinVar.close()
    print("Results are in 'results-ClinVar-long.csv'\n")

print("\n\n\t*** NCBI Search successful ***\n\n")
