import gzip

def clean_fastq(file_name, output_file_name):
    removed_reads = 0  # Tracks the number of removed reads
    total_reads = 0    # Tracks the total number of reads

    with gzip.open(file_name, 'rt') as input_file, gzip.open(output_file_name, 'wt') as output_file:
        while True:
            # Read the four lines that make up a FASTQ record
            read_name = input_file.readline().strip()
            sequence = input_file.readline().strip()
            plus_line = input_file.readline().strip()
            quality_scores = input_file.readline().strip()

            if not read_name or not sequence or not plus_line or not quality_scores:
                break  # End of file or incomplete read
            
            total_reads += 1

            # Check if the read is problematic (e.g., short, invalid characters, incomplete)
            if '===' in sequence or len(sequence) < 20:
                print(f"Removing problematic read: {read_name}")
                print(f"{sequence}")
                removed_reads += 1
                continue

            # Ensure all characters in sequence are valid (ACGTN)
            invalid_chars = set(sequence) - set("ACGTN")
            if invalid_chars:
                print(f"Removing read with invalid characters: {read_name} (Invalid: {invalid_chars})")
                removed_reads += 1
                continue

            # Write valid reads to output file
            output_file.write(read_name + '\n')
            output_file.write(sequence + '\n')
            output_file.write(plus_line + '\n')
            output_file.write(quality_scores + '\n')

    print(f"Total reads processed: {total_reads}")
    print(f"Total problematic reads removed: {removed_reads}")


# Usage
clean_fastq('MUT-CR-1_R1.fastq.gz', 'MUT-CR-1_R1-clean.fastq.gz')
clean_fastq('MUT-CR-1_R2.fastq.gz', 'MUT-CR-1_R2-clean.fastq.gz')