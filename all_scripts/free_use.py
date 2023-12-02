import re
import itertools
import subprocess

import numpy
import pandas as pd
import Bio
from Bio import SeqIO

id_list =[]
with open('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_temp/trichome_ids.txt','r') as id_file:
    for line in id_file:
        id_list.append(line)
tab = pd.read_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/sinapis_alba_3_promoter_analysis_with_TAIR_homologues.txt')
trichome_tab = tab[tab['tair_id'].isin(id_list)]
trichome_tab.to_csv('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/promoter_analysis/sinapis_alba_3_promoter_analysis_with_TAIR_homologues_trichome_genes.txt',sep='\t',quoting=None)
#blast = pd.read_table('/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/Phytozome/PhytozomeV13/Salba/v3.1/annotation/BLASTs/Salba_584_v3.1.cds.fa_all_ids_BLAST_blast_to_TAIR10_cdna_20101214_updated_dedup.txt')



# fa = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome_all_swissprot_gene_hits.fa'
# full_out = '/Volumes/sesame/joerecovery/Project_folder/sinapis_assembly_shenanigans/sinapis_all_rnaseq_reads/maingenome/swissprot_blastx/all_rna_seq_merged_maingenome_tair_swissprot_gene_hits.fa'

# with open(fa,'r') as f, open(full_out,'w+') as out:
#     for seq in SeqIO.parse(f,'fasta'):
#         if 'ARATH' in seq.description:
#             out.write('>' + str(seq.description) + '\n' + str(seq.seq) +'\n')


# def fasta_condenser(fasta, tair=0):
#     '''Preferably use this on the fasta file that has less info/annotation
#     with it - usually the query fasta'''
#     namelist = []
#     seqlist = []
#     if tair==0:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.id)
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['ID', 'Seq'])
#     else:
#         with open (fasta,'r') as fa:
#             for seq in SeqIO.parse(fa,'fasta'):
#                 namelist.append(seq.id[0:9])
#                 seqlist.append(str(seq.seq))
#         df = pd.DataFrame(list(zip(namelist, seqlist)),columns =['TAIR_ID', 'TAIR_Seq'])
#     return df

# def find_nth(haystack, needle, n):
#     start = haystack.find(needle)
#     while start >= 0 and n > 1:
#         start = haystack.find(needle, start+len(needle))
#         n -= 1
#     return start


# id_list = []
# with open('/Users/josephbeegan/Downloads/uniprot-compressed_true_download_true_format_fasta_includeIsoform_tr-2022.11.23-11.40.49.48.fasta') as fasta:
#     for seq in SeqIO.parse(fasta,'fasta'):
#         id_list.append(seq.description)


# real_ids = []
# for i in range(len(id_list)):
#    if 'GN=At' in id_list[i]:
#         real_ids.append(id_list[i])


# swissprot_ids =[]
# tair_ids = []
# for i in real_ids:
#     sp_id_s = find_nth(i, '|', 1)
#     sp_id_f = find_nth(i, '|', 2)
#     tair_id = i.find('GN=At')
#     swissprot_ids.append(i[sp_id_s:sp_id_s + 15])
#     tair_ids.append(i[tair_id+3:tair_id +12].upper())


# for i in range(len(swissprot_ids)):
#     swissprot_ids[i] = swissprot_ids[i][1:10]
# print(swissprot_ids[0:15])
# dict = {'tair': tair_ids, 'swissprot': swissprot_ids}
# # df = pd.DataFrame(dict)
# ids = []
# with open('/Volumes/sesame/joerecovery/genomes/tair_id_list.txt','r') as tair:
#     for line in tair:
#         ids.append(line.strip())
# print(ids)

# with open('/Volumes/sesame/joerecovery/genomes/TAIR_SWISSPROT_IDS.txt','r') as text, open('/Volumes/sesame/joerecovery/genomes/TAIR_SWISSPROT_IDS_fixed.txt','w+') as out:
#     for line in text:
#         if line[0:9] in ids:
#             out.write(line)



# '''
# #Read in tables as pandas dataframes
# blast_output = pd.read_table('MYFILEHERE')
# #Read in the tair to swissprot IDs
# tair_to_swissprot_table = pd.read_table('MYFILEHERE')




# new_blastoutput =  blast_output[blast_output['SWISSPROT_ID'].str.isin('TAIR_ID_LIST')]

# '''
# #print(aug_gtf)
# #print(p_gtf)

# '''



# '''
# # real_list = []
# # for i in genelist:
# #     real_list.append(i[:-2])
# # print(real_list)
# # list = (gtf[gtf.columns[8]].str[9:20])
# # id_list = []
# # for i in list:
# #     id_list.append(i.split('"')[0])
# # gtf['plainids'] = id_list
# # new_gtf = gtf[gtf['plainids'].isin(real_list)]
# # print(new_gtf)
# # print(gtf)
# # '''



# # # def rearrange_hitdata(table):

# # #     outputname = table.split('/')[-1][:-4]
# # #     outputname = outputname + '_rearranged.txt'
# # #     outputpath = table2.split('/')[:-1]
# # #     outputpath = '/'.join(outputpath)
# # #     outputpath = outputpath + '/' + outputname
# # #     table = pd.read_table(table)
# # #     table = table[['Query','Short name']]
# # #     print(outputpath)

# # #     list_for_empty_df = ['Query','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25']
# # #     empty = pd.DataFrame(columns=['Query', '', 'Action'])

# # #     setname = ''
# # #     list_of_ID_names = []
# # #     for index, row in table.iterrows():
# # #         if row['Query'] != setname:
# # #             setname = str(row['Query'])
# # #             list_of_ID_names.append(setname)

# # #     ID_list = table['Query'].tolist()

# # #     max_number = 0
# # #     for i in list_of_ID_names:
# # #         if ID_list.count(i) > max_number:
# # #             max_number= ID_list.count(i)
# # #     list_for_empty_df_new = list_for_empty_df[:max_number +1]
# # #     print(list_for_empty_df_new)
# # #     emptydf = pd.DataFrame(columns=list_for_empty_df_new)
# # #     print(emptydf)


# # #     biglist = []
# # #     for title in list_of_ID_names:
# # #         list = []
# # #         for index, row in table.iterrows():
# # #             if row['Query'] == title:
# # #                 list.append(row['Short name'])
# # #         list.insert(0,title)
# # #         print(str(len(list)))
# # #         biglist.append(list)


# # #     df = pd.DataFrame(biglist,columns = list_for_empty_df_new)
# # #     df.to_csv(outputpath, sep = '\t',index = False)

# # # rearrange_hitdata('/Volumes/seagatedrive/downloads/TRY_hitdata(1).txt')





# # # contigs_with_several_hits_order = open('/Volumes/seagatedrive/Project_folder/sinapis_assembly_shenanigans/blast_script_runs/SOAP_GSS_reads/contigs_with_several_hits_ordered.fasta','w+')
# # # handle = open(seqfile, "rU")
# # # l = SeqIO.parse(handle, "fasta")
# # # sortedList = [f for f in sorted(l, key=lambda x : x.id)]
# # # for s in sortedList:
# # #    contigs_with_several_hits_order.write(s.description + '\n')
# # #    contigs_with_several_hits_order.write(str(s.seq) + '\n')



# # # for i in range(5000000000):
# # #     print('Hi Aisling')
# # # TAIR_gff = pd.read_csv('/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/TAIR10_GFF3_genes.txt', sep='\t')
# # # TAIR_gff = TAIR_gff.drop(['UTR','COPY','TAIR10'],axis =1)
# # # print(TAIR_gff.head())
# # # TAIR_gff = TAIR_gff[TAIR_gff['TYPE']=='protein']
# # # print(TAIR_gff.head())
# # # TAIR_gff['length'] = TAIR_gff['END'] - TAIR_gff['START']
# # # print(TAIR_gff['length'].mean())
# # # print(TAIR_gff['ID2'])
# # # if TAIR_gff['ID2'].str.contains('AT1G01010').any():
# # #     print('wooop')
# # # if 'Bra041148' in 'gene-Bra041148':
# # #     print('yest')
# # # new_blast_SOAP_BO = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/SOAP_ragtag_sinapis_scaffolds_augustus_mrna_BLAST_TAIR10_xml'

# # # yammyg = open('/Volumes/seagatedrive/cLoops/cLoops_env.yaml','r')
# # # conda_installation = open('/Users/josephbeegan/Desktop/cloopspackages.txt','w+')
# # # for line in yammyg:
# # #     if line.startswith(' '):
# # #         line = line[3:]
# # #         full_boi = line.split('=')[0] + line.split('=')[1]
# # #         conda_installation.write(full_boi + ' ')

# # # old_SOAP_blast = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/soapdenovo_JOIGenomics_augustus_mrna_blast.xml'
# # # SOAP_blast_full_assembly = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/SOAP_ragtag_sinapis_scaffolds_full_assembly_BLAST_TAIR10_xml'
# # # soap_augustus_genes = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/SOAP_ragtag_sinapis_scaffolds_augustus.mrna'
# # # tair_genes = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/TAIR10_cdna_20101214_updated'
# # # SOAPvSPADes_blast = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/group_meeting_prep/SPAdes_augustus_to_SOAP_ragtag_augutus_blast_xml'
# # # with open(tair_genes,'r') as handle:
# # #     length_count = 0
# # #     total_count = 0
# # #     for seq in SeqIO.parse(handle,'fasta'):
# # #         print(len(seq.seq))
# # #         length_count = length_count + len(seq.seq)
# # #         total_count +=1
# # #     print(length_count/total_count)
# # # with open(old_SOAP_blast) as soapy:
# # #     soapy_record = NCBIXML.parse(soapy)
# # #     for record in soapy_record:
# # #         for description in record.descriptions:
# # #             for alignment in record.alignments:
# # #                 for hsp in alignment.hsps:
# # #                     if 'AT4G18960' in alignment.title:
# # #                         if description.e < 0:
# # #                             print(record.query)
# # # 135670229
# # # 554977060
# # # SOAP_no_ragtag = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAP_denovo_genomic_assembly'
# # # SOAP_ragtag = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAP_sinapis_assembly_ragtag/SOAP_ragtag_output_no_C_option/SOAP_ragtag_sinapis_scaffolds.fa'
# # # SPAdes_no_ragtag = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/scaffolds.fasta'
# # # SPAdes_ragtag = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/SPAdes_sinapis_assembly_ragtag/ragtag_output/SPAdes_genomic_sinapis_ragtag_scaffolds.fa'
# # # SOAP_59mer = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/SOAPdenvo59mer_run_ragtag_no_correction/ragtag_output/ragtag.scaffolds.fasta'


# # # def seq_count(fasta):
# # #     count = 0
# # #     what_to_print = fasta.split('/')[-1]
# # #     what_to_print = what_to_print[:-3]
# # #     with open(fasta,'r') as handle:
# # #         for seq in SeqIO.parse(handle,'fasta'):
# # #             count += 1
# # #     print(what_to_print+ ' = ' + str(count))
# # # seq_count(SOAP_no_ragtag)
# # # seq_count(SOAP_ragtag)
# # # seq_count(SPAdes_no_ragtag)
# # # seq_count(SPAdes_ragtag)
# # # seq_count(SOAP_59mer)

# # # blast_file = '/Volumes/seagatedrive/emmanuelle_brapa_ids/napus_tair10_cdna_blast_xml'

# # # # snapfasta = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_snap_prediction.txt','r')
# # # # trinity = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/Trinity.fasta.masked','r')
# # # # spades = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdesRNA_assembly.fasta','r')
# # # # contigs = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta','r')
# # # # blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/sinapis_augustus_mrna_tairxml','r')
# # # # trichome_names_tair = ['AT5G53200', 'AT3G27920','AT1G79840','AT2G46410']
# # # # outfile = open('/Users/josephbeegan/Desktop/sinapis_contigs_no_newline.fasta','w+')
# # # #TAIR_genome = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/TAIR10_Chr_all.fasta','r')
# # # # soap_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_augustus.mrna')
# # # # spades_mrna = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/augustus_SPAdes_genomic_assembly/augustu_contigs_Ns.mrna','r')
# # # # spades_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_snap_blast_xml','r')
# # # # spades_snap_genes = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdes_snap.fasta','r')
# # # # soap_snap_blast = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soapdenovo_JOIGenomics_snap_blast_xml','r')
# # # SOAP_ragtag_blast = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_soap_augustus_mrna_tairblast_xml_with_IDs'
# # # TAIR_gff = pd.read_csv('/Volumes/seagatedrive/Work_folder/microarray_SUP/temp_SUP/TAIR10_GFF3_genes.txt', sep='\t')
# # # SPAdes_blast = '/Users/josephbeegan/Desktop/Work/Sinapis_alba/sinapis_augustus_mrna_tairxml'
# # # SPAdes_ragtag_blast = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_augustus_mrna_tairblast_xml_with_IDs'
# # # ZHANG_ragtag_blast = '/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_zhang_augustus_mrna_tairblast_xml_with_IDs'
# # # new_chromosomes= ['Chrom', 'type','genetype','start','end','number','strand','stop','ID']
# # # brassica_gff = pd.read_table('/Volumes/seagatedrive/emmanuelle_brapa_ids/Brapa_gene_v1.5.gff',names = new_chromosomes)
# # # brassica_genome = '/Volumes/seagatedrive/emmanuelle_brapa_ids/newest_fasta_short_ids'


# # # brassica_gff.ID=brassica_gff.ID.str[-9:]
# # # brassica_gff['genetype'] = brassica_gff['genetype'].astype(str) + '-' + brassica_gff['ID'].astype(str)
# # # print(brassica_gff.head())
# # # brassica_gff.to_csv('/Volumes/seagatedrive/edited_brassica_v1.5.gff', sep='\t', index = False)
# # # brassica_gff['test_column'] = brassica_gff['start'].astype(str) + '-' + brassica_gff['end'].astype(str)
# # # brassica_id_list = brassica_gff['ID'].tolist()
# # # new_brassica_gff = brassica_gff[brassica_gff['genetype'] == 'gene']
# # # print(new_brassica_gff.head())
# # # #brassica_gff.ID = brassica_gff.ID.str.split('=')[1]
# # # gene_locations = new_brassica_gff['test_column'].tolist()
# # # gene_ids = new_brassica_gff['ID'].tolist()
# # # print(gene_locations)
# # # print(len(gene_locations))
# # # print(gene_ids)
# # # print(len(gene_ids))
# # # napus_genes = open('/Volumes/seagatedrive/emmanuelle_brapa_ids/wanted_napus_ids_s003.txt','r')
# # # napus_id_list = []
# # # for line in napus_genes:
# # #     napus_id_list.append(line)
# # # print('length' + str(len(line)))
# # # print(napus_id_list)
# # # new_napus_id_list = []
# # # for name in napus_id_list:
# # #     name2 = name[0:9]
# # #     new_napus_id_list.append(name2)


# # # data = []
# # # with open(blast_file) as blast_file2:
# # #     blast_file_record = NCBIXML.parse(blast_file2)
# # #     for nap in new_napus_id_list:
# # #         for record in blast_file_record:
# # #             for description in record.descriptions:
# # #                 for alignment in record.alignments:
# # #                     for hsp in alignment.hsps:
# # #                         # new_napus_summary_file.write('Query ID: ' + record.query + '\n')
# # #                         # new_napus_summary_file.write('Match ID: ' + alignment.title.split(' ')[1]+ '\n')
# # #                         # new_napus_summary_file.write('Query Sequence: '+hsp.query+ '\n')
# # #                         # new_napus_summary_file.write('Score: ' + str(hsp.score)+ '\n')
# # #                         # new_napus_summary_file.write('E-value: ' + str(description.e)+ '\n'+ '\n')
# # #                         data.append([record.query,alignment.title.split(' ')[1],str(description.e)])
# # # s003_df = pd.DataFrame(data, columns=['Brapa_ID', 'TAIR_ID', 'e.value'])

# # #     s003_df = s003_df.drop_duplicates()
# # # print(s003_df.tail())
# # # s003_df.to_csv('/Volumes/seagatedrive/emmanuelle_brapa_ids/brapa_AGI_s002_column_file.txt',index = False, sep = '\t')

# # # for bra in brassica_id_list:
# # #     for nap in new_napus_id_list:
# # #         seen_pair = []
# # #         if bra == nap:
# # #             if bra not in seen_pair:
# # #                 seen_pair.append(bra)
# # #         else:
# # #             continue
# # # print(len(seen_pair))
# # # print(seen_pair)
# # # print(napus_id_list)
# # # print(brassica_id_list)
# # # print(new_napus_id_list)
# # # seqids = []
# # # napus_wanted_fasta = open('/Volumes/seagatedrive/s003_napus_fasta','w+')
# # # with open('/Volumes/seagatedrive/emmanuelle_brapa_ids/newest_fasta_short_ids') as fasta:
# # #     for seq in SeqIO.parse(fasta,'fasta'):
# # #         seqids.append(seq.id)
# # #         for i in range(len(new_napus_id_list)):
# # #             if seq.id == new_napus_id_list[i]:
# # #                 print('yes')
# # #                 napus_wanted_fasta.write('>' + str(seq.id) +'\n'+str(seq.seq)+'\n')
# # # print(seqids)
# # # print(new_napus_id_list)
# # # print(len(set(seqids) & set(new_napus_id_list)))
# # # # outfasta = open('/Volumes/seagatedrive/emmanuelle_brapa_ids/fasta_for_blast','w+')
# # # # with open(brassica_genome) as handle:
# # # #     for seq in SeqIO.parse(handle,'fasta'):
# # # #         for i in range(len((gene_locations))):
# # # #             if gene_locations[i] in seq.id:
# # # #                 print('woo')
# # # #                 outfasta.write(">"+str(gene_ids[i])+'\n' + str(seq.seq)+'\n')

# # # print(brassica_gff.head())


# # # #non_ragtag_spades_file = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/non_ragtag_spades_blast_summary','w+')
# # # #soap_summary_file = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_soap_augustus_mrna_tairblast_quick_summary.txt','w')
# # # #spades_summary_file = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_spades_augustus_mrna_tairblast_quick_summary.txt','w')
# # # zhang_summary_file = open('/Volumes/seagatedrive/sinapis_assembly_shenanigans/ragtag_zhang_augustus_mrna_tairblast_quick_summary.txt','w')
# # # new_blast_SOAP_BO = '/Volumes/seagatedrive/group_meeting_prep/SOAP_ragtag_sinapis_scaffolds_augustus_mrna_BLAST_TAIR10_xml'

# # # with open(new_blast_SOAP_BO) as soapy:
# # #     soapy_record = NCBIXML.parse(soapy)
# # #     for record in soapy_record:
# # #         for description in record.descriptions:
# # #             for alignment in record.alignments:
# # #                 for hsp in alignment.hsps:
# # #                     zhang_summary_file.write('Title: ' + alignment.title + '\n')
# # #                     zhang_summary_file.write('Query: ' + hsp.query + '\n')
# # #                     zhang_summary_file.write('Score: ' + str(description.score) + '\n')
# # #                     zhang_summary_file.write('identity: ' + str(hsp.identities) + '\n')
# # #                     zhang_summary_file.write('E-score: ' + str(description.e) + '\n')
# # #                     zhang_summary_file.write('Bitscore: ' + str(hsp.bits) + '\n' + '\n')



# # # for seq in SeqIO.parse(TAIR_genome,'fasta'):
# # #     print(str(seq.seq)[12562141:12562163])
# # #     primer = 'TCCGACAAGATCCTTTGAGCAC'
# # #     primer2 = Seq('ATAAATATGGATTGCACGACGAGTC')
# # #     if primer in seq.seq:
# # #         print(primer + ' Found in ' + seq.id)
# # #         print(re.search(primer, str(seq.seq)).start())
# # #     else:
# # #         print(primer + ' Not in ' + seq.id)
# # #     if str(primer2.reverse_complement()) in seq.seq:
# # #         print( primer2 + ' Found in ' + seq.id)
# # #     else:
# # #         print(primer2 + ' Not in ' + seq.id)
# # # line = ''
# # # bamming_command = ['echo', 'fart', line]
# # # file = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/queries_file.txt','r')
# # # for line2 in file:
# # #     line = line2
# # #     subprocess.Popen(bamming_command, stdout=subprocess.PIPE, bufsize=-1)

# # # f2 = pysam.AlignmentFile('/Volumes/seagatedrive/sinapis_assembly_shenanigans/SRR7224588_R1_10000_reads.sam','r')
# # # list = []
# # # full_count = 0
# # # unmapped_count = 0
# # # for read in f2:
# # #     full_count +=1
# # #     if read.mapping_quality > 0:
# # #         list.append(read)
# # #     else:
# # #         unmapped_count += 1
# # # print('There are {} mapped reads and {} unmapped reads out of {}'.format(len(list),str(unmapped_count),str(full_count)))

# # # def levenshtein_ratio_and_distance(s, t, ratio_calc = False):
# # #     """ levenshtein_ratio_and_distance:
# # #         Calculates levenshtein distance between two strings.
# # #         If ratio_calc = True, the function computes the
# # #         levenshtein distance ratio of similarity between two strings
# # #         For all i and j, distance[i,j] will contain the Levenshtein
# # #         distance between the first i characters of s and the
# # #         first j characters of t
# # #     """
# # #     # Initialize matrix of zeros
# # #     rows = len(s)+1
# # #     cols = len(t)+1
# # #     distance = np.zeros((rows,cols),dtype = int)

# # #     # Populate matrix of zeros with the indeces of each character of both strings
# # #     for i in range(1, rows):
# # #         for k in range(1,cols):
# # #             distance[i][0] = i
# # #             distance[0][k] = k

# # #     # Iterate over the matrix to compute the cost of deletions,insertions and/or substitutions
# # #     for col in range(1, cols):
# # #         for row in range(1, rows):
# # #             if s[row-1] == t[col-1]:
# # #                 cost = 0 # If the characters are the same in the two strings in a given position [i,j] then the cost is 0
# # #             else:
# # #                 # In order to align the results with those of the Python Levenshtein package, if we choose to calculate the ratio
# # #                 # the cost of a substitution is 2. If we calculate just distance, then the cost of a substitution is 1.
# # #                 if ratio_calc == True:
# # #                     cost = 2
# # #                 else:
# # #                     cost = 1
# # #             distance[row][col] = min(distance[row-1][col] + 1,      # Cost of deletions
# # #                                  distance[row][col-1] + 1,          # Cost of insertions
# # #                                  distance[row-1][col-1] + cost)     # Cost of substitutions
# # #     if ratio_calc == True:
# # #         # Computation of the Levenshtein Distance Ratio
# # #         Ratio = ((len(s)+len(t)) - distance[row][col]) / (len(s)+len(t))
# # #         return Ratio
# # #     else:
# # #         # print(distance) # Uncomment if you want to see the matrix showing how the algorithm computes the cost of deletions,
# # #         # insertions and/or substitutions
# # #         # This is the minimum number of edits needed to convert string a to string b
# # #         return "The strings are {} edits away".format(distance[row][col])

# # # def blast_to_dataframe(BLAST_record):
# # #     titles = []
# # #     starts = []
# # #     lengths = []
# # #     matches = []
# # #     queries = []
# # #     for blast_record in BLAST_record:
# # #         for description in blast_record.descriptions:
# # #             for alignment in blast_record.alignments:
# # #                 for hsp in alignment.hsps:
# # #                     '''From here it is only taking the trichome gene_IDs'''
# # #                     for name in trichome_names_tair:
# # #                         if name in alignment.title:
# # #                             titles.append(alignment.title.split(' ')[1])
# # #                             starts.append(hsp.query_start)
# # #                             lengths.append(alignment.length)
# # #                             queries.append(hsp.query.split('-'))



# # #     pre_dict = {'titles':titles,'lengths':lengths,'starts':starts,'queries':queries}
# # #     datafram = pd.DataFrame(pre_dict)
# # #     return(datafram)

# # # record = NCBIXML.parse(blast)
# # # soaprecord = NCBIXML.parse(soap_blast)
# # # soap_snap_record = NCBIXML.parse(soap_snap_blast)
# # # spades_snap_record = NCBIXML.parse(spades_snap_blast)




# # # soap_augustus = pd.read_table('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_augustus_trichome_genes_full_df')
# # # spades_augustus = pd.read_table('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_augustus_trichome_genes_full_df')
# # # soap_snap = pd.read_table('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_snap_trichome_genes_full_df')



# # # augustus_merge = soap_augustus.merge(spades_augustus,right_on = 'TAIR_titles',left_on='TAIR_titles')
# # # print(augustus_merge.head())
# # # print(augustus_merge.columns)
# # # print(augustus_merge)
# # # augustus_merge.to_csv('/Users/josephbeegan/Desktop/soap_augustus_spades_augutus_merged_table',index = False, sep = '\t')
# # # soap_pre = ''
# # # spades_pre = ''
# # # print(augustus_merge[['TAIR_titles','Augustus_titles_x','Augustus_titles_y']])
# # # spades_count = 0


# # # spades_augustus['predicted_genes'] = ''
# # # for seq in SeqIO.parse(spades_mrna,'fasta'):
# # #     for i in range(len(spades_augustus)):
# # #         if spades_augustus['Augustus_titles'][i] in seq.id:
# # #             spades_augustus['predicted_genes'][i] = str(seq.seq)
# # # print(spades_augustus[['Augustus_titles','predicted_genes']])


# # # soap_augustus['predicted_genes'] = ''
# # # for seq in SeqIO.parse(soap_mrna,'fasta'):
# # #     for i in range(len(soap_augustus)):
# # #         if soap_augustus['Augustus_titles'][i] in seq.id:
# # #             soap_augustus['predicted_genes'][i] = str(seq.seq)
# # # print(soap_augustus[['Augustus_titles','predicted_genes']])
# # # #spades_augustus_tair_soap_match = spades_augustus['TAIR_titles'] in soap_augustus['TAIR_titles']
# # # #print(spades_augustus[spades_augustus_tair_soap_match])
# # # trimmed_spades_augustus = spades_augustus[['TAIR_titles','Augustus_titles','predicted_genes']]
# # # trimmed_soap_augustus = soap_augustus[['TAIR_titles','Augustus_titles','predicted_genes']]
# # # print(trimmed_soap_augustus.head())
# # # print(trimmed_spades_augustus.head())
# # # trimmed_augustus_merge = trimmed_soap_augustus.merge(trimmed_spades_augustus,right_on = 'TAIR_titles',left_on='TAIR_titles')
# # # for i in range(len(trimmed_augustus_merge)):
# # #     print(trimmed_augustus_merge['TAIR_titles'][i])
# # #     print(str(levenshtein_ratio_and_distance(trimmed_augustus_merge['predicted_genes_x'][i],trimmed_augustus_merge['predicted_genes_y'][i],ratio_calc = True)))



# # # def compare_gene_predictions(df,soapmrna,spadesmrna,TAIRID):
# # #     spades_pre = ''
# # #     for spadesseq in SeqIO.parse(spadesmrna,'fasta'):
# # #         spades_break = 0
# # #         for i in range(len(df)):
# # #             if TAIRID in df['TAIR_titles'][i]:
# # #                 if spadesseq.id == df['Augustus_titles_y'][i]:
# # #                     spades_pre = spadesseq.seq
# # #                     print(spadesseq.id)
# # #             else:
# # #                 spades_break = 1
# # #     soap_pre = ''
# # #     for soapseq in SeqIO.parse(soapmrna,'fasta'):
# # #         soap_break = 0
# # #         for i in range(len(df)):
# # #             if TAIRID in df['TAIR_titles'][i]:
# # #                 if soapseq.id == df['Augustus_titles_x'][i]:
# # #                     soap_pre = soapseq.seq
# # #                     print(soapseq.id)
# # #             else:
# # #                 soap_break = 1
# # #     if len(spades_pre) and len(soap_pre) > 0:
# # #         print('Comparing ' + str(TAIRID) + ': ' + str(levenshtein_ratio_and_distance(soap_pre,spades_pre,ratio_calc = True)))
# # #     if soap_break or spades_break > 0:
# # #         print('Only one of the predictions had a sequence matching ' + str(TAIRID))
# # # for i in range(len(augustus_merge)):
# # #     TAIR = str(augustus_merge['TAIR_titles'][i])
# # #     print(TAIR)
# # #     compare_gene_predictions(augustus_merge,soap_mrna,spades_mrna,TAIR)
# # # #for i in range(len(spades_augustus)):
# # # #    if spades_augustus['TAIR_titles'] == soap_augustus['TAIR_titles'].values:
# # # #        print('wwoppppopo')

# # # for soap in soap_augustus.values:
# # #     print(soap)
# # #             for soap_seq in (soap_augustus['queries'][i].split(',')):
# # #                 for spades_seq in (spades_augustus['queries'][j].split(',')):
# # #                     print('Comparing ' + str(soap_augustus['SOAP_titles'][j]) + ' to ' + str(spades_augustus['SPAdes_titles'][i]) + ':' +levenshtein_ratio_and_distance(soap_seq,spades_seq))

# # # for seq in SeqIO.parse(contigs,'fasta'):
# # #     for i in range(len(spades_augustus)):
# # #         for que in spades_augustus['queries'][i]:
# # #             if len(que) > 14:
# # #                 if que.upper() in seq.seq.upper():
# # #                     print(que)
# # #                     test_dict.setdefault(seq.id + '_' + dataframe['TAIR_titles'][i],[])
# # #                     test_dict[seq.id].append(que)



# # # #soap_snap_df = blast_to_dataframe(soap_snap_record)
# # # #soap_augustus_df = blast_to_dataframe(soaprecord)
# # # #spades_snap_df = blast_to_dataframe(spades_snap_record)
# # # #spades_augutus_df = blast_to_dataframe(record)


# # # print(make_query_dict(soap_augustus))




# # # soap_augustus = pd.read_table('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_augustus_trichome_genes_full_df')
# # # spades_augustus = pd.read_table('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_augustus_trichome_genes_full_df')
# # # probabilities = {}
# # # for j in range(len(soap_augustus)):
# # #     probabilities[[str(soap_augustus['TAIR_titles'][j]]] = str(soap_augustus['queries'][j])
# # # print(probabilities)
# # # for line in str(soap_augustus['queries']).split(','):
# # #     print(line)
# # # print(dict)




















# # # def add_contig_sequences_to_df(df,contigs,contig_title):
# # #     df['sequence'] = ''
# # #     for seq in SeqIO.parse(contigs,'fasta'):
# # #         for i in range(len(df)):
# # #             if seq.id == df[contig_title][i]:
# # #                 df['sequence'] = str(seq.seq)
# # #     return(df)
# # # #for spades in SeqIO.parse(spades_contigs,'fasta'):
# # # #    for i in range(len(spades_augustus_table)):
# # # #        if spades.id == spades_augustus_table['SPAdes_title'][i]:
# # # #            spades_augustus_table['sequence'][i] = str(spades.seq)
# # # spades_contigs = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta','r')
# # # soap_contigs = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SOAP_denovo_genomic_assembly','r')
# # # spades_augustus = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SPAdes_augustus_genomic_assembly_trichome_matches','r')
# # # soap_augustus = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SOAPdenovo_augustus_genomic_assembly_trichome_matches','r')
# # # soap_snap = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/SOAPdenovo_snap_genomic_assembly_trichome_matches','r')
# # # spades_contigs_IDS = []
# # # spades_contigs_seqs = []
# # # for seq in SeqIO.parse(spades_contigs,'fasta'):
# # #     spades_contigs_IDS.append(seq.id)
# # #     spades_contigs_seqs.append(str(seq.seq))
# # # spades_contigs_df = pd.DataFrame(list(zip(spades_contigs_IDS,spades_contigs_seqs)),columns = ['seq.id','seq.seq'])
# # # soap_contigs_IDS = []
# # # soap_contigs_seqs = []
# # # for seq in SeqIO.parse(soap_contigs,'fasta'):
# # #     soap_contigs_IDS.append(seq.id)
# # #     soap_contigs_seqs.append(str(seq.seq))
# # # soap_contigs_df = pd.DataFrame(list(zip(soap_contigs_IDS,soap_contigs_seqs)),columns = ['seq.id','seq.seq'])
# # # spades_augustus_table = pd.read_table(spades_augustus)
# # # soap_augustus_table = pd.read_table(soap_augustus)
# # # soap_snap_table = pd.read_table(soap_snap)


# # # print(spades_augustus_table.head())
# # # print(spades_augustus_table['TAIR_titles'])
# # # merged_spades_augutus = spades_augustus_table.merge(spades_contigs_df, left_on='SPAdes_title', right_on='seq.id')
# # # print(merged_spades_augutus.columns)
# # # print(merged_spades_augutus.head())
# # # print(merged_spades_augutus['SPAdes_title'])
# # # print(merged_spades_augutus['seq.id'])
# # # print(merged_spades_augutus['seq.seq'])
# # # merged_spades_augutus.drop(['seq.id'],axis=1)
# # # merged_spades_augutus.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_augustus_trichome_genes_full_df',sep='\t',index = False)

# # # merged_soap_augustus = soap_augustus_table.merge(soap_contigs_df,left_on='SOAP_titles',right_on='seq.id')
# # # merged_soap_augustus.drop(['seq.id'],axis = 1)
# # # merged_soap_augustus.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_augustus_trichome_genes_full_df',index=False,sep='\t')
# # # print(merged_soap_augustus.head())
# # # merged_soap_snap = soap_snap_table.merge(soap_contigs_df,left_on='SOAP_titles',right_on='seq.id')
# # # merged_soap_snap.drop(['seq.id'],axis = 1)
# # # merged_soap_snap.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_snap_trichome_genes_full_df',index=False,sep='\t')
# # # print(merged_soap_snap.head())
# # # #for spades in SeqIO.parse(spades_contigs,'fasta'):
# # # #    for i in range(len(spades_augustus_table)):
# # # #        if spades.id == spades_augustus_table['SPAdes_title'][i]:
# # # #            spades_augustus_table['sequence'][i] = str(spades.seq)
# # # add_contig_sequences_to_df(spades_augustus_table,spades_contigs,'SPAdes_title')
# # # add_contig_sequences_to_df(soap_augustus_table,soap_contigs,'SOAP_titles')
# # # add_contig_sequences_to_df(soap_snap_table,soap_contigs,'SOAP_titles')
# # # print(spades_augustus_table.head())



# # # spades_augustus_outfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/spades_augustus_contigs_with_matches','w+')
# # # soap_augustus_outfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_augustus_contigs_with_matches','w+')
# # # soap_snap_outfile = open('/Users/josephbeegan/Desktop/Work/Sinapis_alba/soap_snap_contigs_with_matches','w+')
# # # for i in range(len(spades_augustus_table)):
# # #     spades_augustus_outfile.write(spades_augustus_table['TAIR_titles'][i] + '_' + spades_augustus_table['SPAdes_title'][i])
# # #     spades_augustus_outfile.write(spades_augustus_table['sequences'][i])
# # # expected_result = pd.merge(spades_augustus_table, soap_augustus_table, on = 'TAIR_titles', how = 'left')
# # # print(expected_result.head())
# # # expected_result.to_csv('/Users/josephbeegan/Desktop/Work/Sinapis_alba/matches_between_spades_and_soap_augustus', sep = '\t',index=False)

# # # for i in range(len(spades_augustus_table)):
# # #     if spades_augustus_table['titles'][i]
