`import tempfile
import logging
import os
import pysam
import subprocess
fart = '/Volumes/LaCie/READS/joe/Hi-C/raw_reads_SRR7224588/raw_reads1/SRR7224588_R1_10000_reads.fastq'
print(fart.split())
'''Map raw HiC reads iteratively with bowtie2. http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml Iterative mapping accounts for the modification of fragments’ sequences due to ligation.'''

'''This function is very straightforward, with checking if the mapping quality of the read in the outsamfile is >0 we can see if it was mapped at all,
this will eventually be put through another function to extract the unmapped read-ids which will then be elongated and mapped again'''
MIN_MAPQ = 1
def readIsUnmapped(read):
    if (read.mapping_quality < MIN_MAPQ):
        return True

    # Skip non-uniquely aligned.
    for tag in read.tags:
        if tag[0] == 'XS':
            return True
    return False

def sleep():
    """sleep for a second, run garbage collector, sleep again.
    Sleep is split in small pieces to allow some callbacks to
    possibly terminate in between (I don't know if it makes sense, but
    it definitely does not hurt)"""
    for _ in range(3):
        time.sleep(0.1)
    gc.collect()
    for _ in range(3):
        time.sleep(0.1)

'''The algorithm of iterative correction:

    Truncate the sequences to the first N = min_seq_len base pairs, starting at the seq_start position.
    Map the sequences using bowtie2.
    Store the uniquely mapped sequences in a SAM file at out_sam_path.
    Go to the step 1, increase the truncation length N by len_step base pairs, and map the non-mapped and non-uniquely mapped sequences, ... Stop when the 3’ end of the truncated sequence reaches the seq_end position.
'''
def _filter_fastq(ids, inStream, out_fastq, in_filename="none"):  # @UnusedVariable
    '''Filter FASTQ sequences by their IDs.

    Read entries from **in_fastq** and store in **out_fastq** only those
    the whose ID are in **ids**.
    '''
    writingProcess = open(out_fastq.'w+')

    num_filtered = 0
    num_total = 0
    while True:

        line = inStream.readline()
        try:
            assert six.indexbytes(line,0) == 64 # "@"
        except AssertionError:
            print('Not fastq')
            raise
        except IndexError:
            break


        # raise Exception('{0} does not comply with the FASTQ standards.'.format(in_filename))

        fastq_entry = (line, inStream.readline(),
                       inStream.readline(), inStream.readline())
        read_id = line.split()[0][1:]
        if read_id in ids:
            writingProcess.stdin.writelines(fastq_entry)
            num_filtered += 1
        num_total += 1


    sleep()
    writingProcess.communicate()

    if writingProcess.returncode != 0:
        raise RuntimeError("Writing process return code {0}".format(writingProcess.returncode))
    return num_total, num_filtered
def _filter_unmapped_fastq(in_stream, in_sam, nonunique_fastq, in_filename="none"):
    '''Read raw sequences from **in_fastq** and alignments from
    **in_sam** and save the non-uniquely aligned and unmapped sequences
    to **unique_sam**.
    '''
    #in_sam in this case would be the
    samfile = pysam.Samfile(in_sam)  # @UndefinedVariable
    nonunique_ids = set()
    for read in samfile:
        if readIsUnmapped(read):
            nonunique_ids.add(read.qname.encode())

    num_total, num_filtered = _filter_fastq(
        nonunique_ids, in_stream, nonunique_fastq, in_filename=in_filename)
    sleep()

    return num_total, num_filtered



def iterative_mapping(fastq_path, out_sam_path,
                      min_seq_len, len_step, **kwargs):
    bowtie_path =  '/Users/frankwellmer/miniconda3/envs/bio_p3/bin/bowtie2'
    bowtie_index_path = '/Users/frankwellmer/refs/Bowtie2Index/genome'
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)

    '''This is the general command to use, first need to figure out how to trim the reads to the right length - the bowtie2 options -5 and
    -3 allow for trimming of n number of bases from each read before alignment from the 5' end and the 3' end, respectively '''
    #mapping_command = [
    #            bowtie_path, '-x', bowtie_index_path, '-q', '-',
    #            '-5', str(trim_5), '-3', str(trim_3), '-p', str(nthreads)
    #            ] + bowtie_flags.split()
    seq_start = kwargs.get('seq_start', 0)
    seq_end = kwargs.get('seq_end', None)
    nthreads = kwargs.get('nthreads', 8)
    max_len = kwargs.get("max_len", 9999)
    log.info("Using new argument: max_len = {0}".format(max_len))
    bowtie_flags = kwargs.get('bowtie_flags', '')
    fastq_path = os.path.abspath(os.path.expanduser(fastq_path))
    if not os.path.isfile(fastq_path):
        raise Exception(
            'The fastq file is not found '
            'at the specified path: {0}.'.format(fastq_path))
    print(fastq_path)
    extension = fastq_path.split('.')[-1].lower()
    if extension == 'gz':
        bash_reader = 'gunzip -c'
    if extension == 'fastq':
        bash_reader = 'cat'
    reading_command = bash_reader.split() + [fastq_path, ]

    print(fastq_path)
    print(reading_command)
    reading_command = bash_reader.split() + [fastq_path, ]
    reading_process = subprocess.Popen(reading_command,
                                           stdout=subprocess.PIPE)
    reading_process.stdout.readline()
    raw_seq_len = len(reading_process.stdout.readline().strip())

    temp_dir = os.path.abspath(os.path.expanduser(
            kwargs.get('temp_dir', tempfile.gettempdir())))
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)


    '''This section here will just drop the PHRED scores from the FASTQ files to save space, they shouldnt really be necessary if the
    quality check has already been done'''

    drop_seqs_command = ['awk',
                """{OFS="\\t"; if ($1 ~ !/^@/) { $10="A"; $11="g"; if ($3 ~ /\\*/) $6="*"; else $6="1M"; } print}"""]
    output_is_bam = (out_sam_path.split('.')[-1].lower() == 'bam')
    bamming_command = ['samtools', 'view', '-bS', '-'] if output_is_bam else []

    '''os.system and subprocess are essentially the same thing. So, if i were to do os.system.call it would do the same thing as
    subprocess.call - apparently subprocess is a better way to go about it'''



    '''This is how were are going to go to the next iteration - pysam is a python module used to read sam files'''

    '''This is the part where the trimming of the reads and the alignment actually occurs'''
    local_seq_end = min(raw_seq_len, seq_end) if seq_end else raw_seq_len

    if min_seq_len <= local_seq_end - seq_start:
        trim_5 = seq_start
        trim_3 = raw_seq_len - seq_start - min_seq_len
        if raw_seq_len - trim_3 - trim_5 > max_len:
            trim_5 = raw_seq_len - trim_3 - max_len
        local_out_sam = out_sam_path + '.' + str(min_seq_len)
        mapping_command = [
            bowtie_path, '-x', bowtie_index_path, '-q', '-',
            '-5', str(trim_5), '-3', str(trim_3), '-p', str(nthreads)
            ] + bowtie_flags.split()
        pipeline = []
        try:

            '''If i am understanding this correctly, we are creating a pipeline in the form of a list
            and then taking the last argument from the list each time (which is the next read to map) then,
            we are taking the PHRED score and other unneccessary columns from the read, then mapping it - writing
            to local_out_sam'''
            log.info('Reading command: %s', ' '.join(reading_command))
            pipeline.append(
            subprocess.Popen(reading_command, stdout=subprocess.PIPE, bufsize=-1))
            log.info('Mapping command: %s', ' '.join(mapping_command))
            pipeline.append(
                subprocess.Popen(mapping_command,
                    stdin=pipeline[-1].stdout,
                    stdout=subprocess.PIPE if (bamming_command or drop_seqs_command) else open(local_out_sam, 'w'),
                    bufsize=-1))

            if drop_seqs_command:
                log.info('Output editing command: %s', ' '.join(drop_seqs_command))
                pipeline.append(
                    subprocess.Popen(drop_seqs_command,
                        stdin=pipeline[-1].stdout,
                        stdout=subprocess.PIPE if bamming_command else open(local_out_sam, 'w'),
                        bufsize=-1))

            if bamming_command:
                log.info('Output formatting command: %s', ' '.join(bamming_command))
                pipeline.append(
                    subprocess.Popen(bamming_command,
                        stdin=pipeline[-1].stdout,
                        stdout=open(local_out_sam, 'w'),
                        bufsize=-1))
            pipeline[-1].wait()
            '''I think this next part is just making sure that the dropping and mapping has occurred (or been
            attempted) for all of the reads'''
        finally:
            for process in pipeline:
                if process.poll() is None:
                    process.terminate()

        if (len_step <= 0) or (min_seq_len + len_step > local_seq_end - seq_start):
            if kwargs.get("first_iteration", True) == False:
                print("Deleting previous file", fastq_path)
                os.remove(fastq_path)
            return
        log.info('Save the unique aligments and send the '
                     'non-unique ones to the next iteration')
        reading_process = subprocess.Popen(reading_command,
                                       stdout=subprocess.PIPE,
                                       bufsize=-1)

        unmapped_fastq_path = os.path.join(
            temp_dir, os.path.split(fastq_path)[1] + '.%d' % min_seq_len + ".fastq.gz")

        num_total, num_filtered = _filter_unmapped_fastq(
            reading_process.stdout, local_out_sam, unmapped_fastq_path, in_filename=fastq_path)

        reading_process.communicate()
        sleep()


        log.info(('{0} non-unique reads out of '
                  '{1} are sent the next iteration.').format(num_filtered, num_total))

        if kwargs.get("first_iteration", True) == False:
            print("Deleting previous file", fastq_path)
            os.remove(fastq_path)

        kwargs['first_iteration'] = False
        if commandExists("pigz"):
            kwargs["bash_reader"] = "pigz -dc"
        else:
            kwargs["bash_reader"] = 'gunzip -c'
        min_seq_len=min_seq_len + len_step,

        iterative_mapping(bowtie_path, bowtie_index_path, unmapped_fastq_path,
                          out_sam_path,
                          min_seq_len=min_seq_len + len_step,
                          len_step=len_step, **kwargs)


iterative_mapping('/Volumes/LaCie/READS/joe/Hi-C/raw_reads_SRR7224588/raw_reads1/SRR7224588_R1.fastq', 'out_sam_iterative_mapping_trial', 25, 5)
