
sinapis_protein_tair_matches = {'Tair_ID':ids,'Sinapis_query':query}


sinapis_protein_tair_matchesdf = pd.DataFrame.from_dict(sinapis_protein_tair_matches)



for query in sinapis_protein_tair_matchesdf['Sinapis_query']:
    for seq in SeqIO.parse(sinapis_transcriptome, 'fasta'):
        if query in seq.seq:
            print(sinapis_protein_tair_matchesdf[sinapis_protein_tair_matchesdf['Sequence']==string]['Length'].values)
        print(pairwise2.align.localxx(seq.seq,query))



def fasta_sequences_recover(fasta):
    sequence_list = []
    for seq in SeqIO.parse(fasta,'fasta'):
        seqs = str(seq.seq)
        sequence_list.append(seqs)
    return(sequence_list)


brapa_hairy_genes = "/Users/josephbeegan/Desktop/Work/Sinapis_alba/brassica_gene_sequences.fa"
brapa_genes = fasta_sequences_recover(brapa_hairy_genes)
brapa_query = brapa_genes[0]
sinapis_alba_genes = "/Users/josephbeegan/Desktop/Work/Sinapis_alba/contigs.fasta"
sinapis_alba_genes = fasta_sequences_recover(sinapis_alba_genes)
sinapis_alba_subject = sinapis_alba_genes[0]

alignments = pairwise2.align.globalxx(sinapis_alba_subject, brapa_query)
print(format_alignment(*alignments[0]))

#for a in pairwise2.align.localxx(sina
pis_alba_subject, brapa_query):
#    if len(format_alignment(*a)) > 500:
#        print(format_alignment(*a))
start = 'ATG'
stop = ['TGA','TAG','TAA']
start_codon_positions = []
stop_codon_positions = []
for i in range(0, len(sinapis_alba_subject)-2, 3):
    codon = sinapis_alba_subject[i:i+3]
    if codon == start:
        start_codon_positions.append(i)
    if codon in stop:
        stop_codon_positions.append(i)
potential_orfs = []
for start in start_codon_positions:
    for stop in stop_codon_positions:
        if start < stop:
            potential_orfs.append(sinapis_alba_subject[start:stop+2])

outfasta_potential_orfs = open('/Users/josephbeegan/Documents/sinapis_alba_potential_orfs.fa','w+')
for orf,i in zip(potential_orfs, range(len(potential_orfs))):
    outfasta_potential_orfs.write('>ORF' + str(i) + '\n')
    outfasta_potential_orfs.write(orf + '\n')




def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, reverse_complement(seq))]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1
    answer.sort()
    return answer

for orf in potential_orfs:
    orf_list = find_orfs_with_trans(orf, table, min_pro_len)
for start, end, strand, pro in orf_list:
    print(
        "%s...%s - length %i, strand %i, %i:%i"
        % (pro[:30], pro[-3:], len(pro), strand, start, end)
    )

test = "AAABBBCCCDDDEEEFFFGGGHHHIII"
for i in range(0, len(test)-2, 3):
    print(test[i:i+3])
