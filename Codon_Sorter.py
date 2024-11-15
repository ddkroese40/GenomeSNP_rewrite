import sys
import re

def sort_fasta_by_gene_number(input_file):
    with open(input_file, 'r') as f:
        content = f.read().split('>')[1:]  # Skip the first empty element from the split

    # Initialize a list to store entries
    fasta_entries = []

    # Parse the fasta content into blocks of entries for each gene
    i = 0
    while i < len(content):
        lines = content[i].strip().splitlines()
        header = lines[0]
        sequence = '\n'.join(lines[1:])

        # Extract the gene number
        gene_number_match = re.search(r'gene=Hg_chrom1_TN10gene_(\d+)', header)
        if gene_number_match:
            gene_number = int(gene_number_match.group(1))
            entry_block = [f'>{header}\n{sequence}']

            # Check and add replicate_1 and replicate_2, if they exist
            for j in range(1, 3):  # look ahead for two replicate entries
                if i + j < len(content):
                    replicate_lines = content[i + j].strip().splitlines()
                    if replicate_lines[0].startswith("replicate_"):
                        entry_block.append(f'>{replicate_lines[0]}\n' + '\n'.join(replicate_lines[1:]))
                    else:
                        break

            fasta_entries.append((gene_number, entry_block))
            i += len(entry_block)  # Move index ahead by the number of lines for this gene + replicates
        else:
            i += 1

    # Sort entries by gene number
    fasta_entries.sort(key=lambda x: x[0])

    # Define the output filename
    output_file = input_file.replace('.fasta', '_sorted.fasta')

    # Write sorted entries to the output file
    with open(output_file, 'w') as f:
        for _, entry_block in fasta_entries:
            f.write('\n'.join(entry_block) + '\n')

    print(f"Sorted file created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python sort_fasta.py <input_file.fasta>")
    else:
        input_file = sys.argv[1]
        sort_fasta_by_gene_number(input_file)
