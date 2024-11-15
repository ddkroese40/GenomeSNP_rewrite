import argparse
import os

def parse_cds_positions(file_path):
    cds_positions = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 5 and columns[3].isdigit() and columns[4].isdigit():  # Ensure start and end are present and numeric
                feature_type = columns[2]
                if feature_type == 'CDS':
                    start = int(columns[3]) - 1  # Convert to zero-indexed
                    end = int(columns[4]) - 1
                    cds_positions.append((start, end))
    return cds_positions


def calculate_cds_coverage(transcript_dict, cds_positions):
    cds_index = 0
    
    # Iterate through each transcript entry
    for rna_name, transcript in transcript_dict.items():
        transcript_length = len(transcript["replicate1"])  # Assume replicate1 and replicate2 are the same length
        cumulative_cds_length = 0
        cds_count = 0
        
        # Accumulate CDS lengths until the cumulative length matches the transcript length
        while cds_index < len(cds_positions):
            start, end = cds_positions[cds_index]
            cds_length = end - start + 1
            cumulative_cds_length += cds_length
            cds_count += 1
            cds_index += 1

            # Check if the cumulative length matches the transcript length
            if cumulative_cds_length == transcript_length:
                break
            elif cumulative_cds_length > transcript_length:
                print("ERROR: CDS length longer than transcript", rna_name, "cumulative length", cumulative_cds_length, "transcript length:", transcript_length)
                break

        # Add the number of CDS regions represented by the transcript to transcript_dict
        transcript_dict[rna_name]["cds_count"] = cds_count
        #if you want to view how many proteins per transcript for each
        #print("index:", rna_name, "count:", cds_count)


def RewriteDNA(file_path, source_gene_path, original_gene_path):
    # Parse CDS positions and store them as a list of tuples
    cds_positions = parse_cds_positions(file_path)
    
    # Parse transcript information from the source gene file
    transcript_dict = ParseTranscript(source_gene_path)

    # Calculate and update CDS coverage for each transcript
    calculate_cds_coverage(transcript_dict, cds_positions)

    # Call ReplaceCDS with the list of CDS positions and the source gene file
    ReplaceCDS(transcript_dict, cds_positions, original_gene_path, file_path)


def replace_cds_for_replicate(cds_slice, original_sequence, replacement_sequence):
    """
    Helper function to replace CDS sequences in the original gene for one replicate.
    """
    updated_sequence = list(original_sequence)
    for start, end in cds_slice:
        updated_sequence[start:end + 1] = list(replacement_sequence[:end - start + 1])
        replacement_sequence = replacement_sequence[end - start + 1:]
    return ''.join(updated_sequence)


def write_replicate_file(file_path, replicate_sequence, replicate_name):
    def format_sequence(sequence):
        return '\n'.join(sequence[i:i + 59] for i in range(0, len(sequence), 59))
    
    # Remove the file extension from file_path
    base_file_name = os.path.splitext(file_path)[0]
    output_file = f"{base_file_name}_{replicate_name}_out.fasta"
    
    # Write the formatted sequence to the output file
    with open(output_file, 'w') as file:
        file.write(f">{replicate_name}\n" + format_sequence(replicate_sequence) + "\n")


def ReplaceCDS(transcript_dict, cds_positions, original_gene_path, file_path):
    # Load the original gene sequence as a single string
    with open(original_gene_path, 'r') as file:
        original_gene = ''.join(line.strip() for line in file if not line.startswith('>'))
    
    replicate1_sequence = original_gene
    replicate2_sequence = original_gene
    cds_index = 0

    # Iterate through each transcript entry in transcript_dict
    for rna_name, transcript in transcript_dict.items():
        replicate1_replacement = transcript["replicate1"]
        replicate2_replacement = transcript["replicate2"]
        cds_count = transcript["cds_count"]

        # Get the current slice of CDS positions
        cds_slice = cds_positions[cds_index:cds_index + cds_count]

        # Replace CDS sequences in the original gene for each replicate
        replicate1_sequence = replace_cds_for_replicate(cds_slice, replicate1_sequence, replicate1_replacement)
        replicate2_sequence = replace_cds_for_replicate(cds_slice, replicate2_sequence, replicate2_replacement)

        # Increment the cds_index by cds_count
        cds_index += cds_count

    # Write the modified sequences to output files
    write_replicate_file(file_path, replicate1_sequence, "Replicate1")
    write_replicate_file(file_path, replicate2_sequence, "Replicate2")


def ParseTranscript(source_gene_path):
    transcript_dict = {}
    current_rna = None
    replicate_type = None  # Initialize replicate_type here to avoid UnboundLocalError
    
    with open(source_gene_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if "replicate_1" in line:
                    replicate_type = "replicate1"
                elif "replicate_2" in line:
                    replicate_type = "replicate2"
                else:
                    # New RNA entry
                    current_rna = line.split()[0][1:]  # Remove the ">" and take the RNA name
                    transcript_dict[current_rna] = {"replicate1": "", "replicate2": ""}
                    replicate_type = None  # Reset replicate_type for new RNA
            elif current_rna and replicate_type:
                # Assign the sequence to the correct replicate
                transcript_dict[current_rna][replicate_type] = line
                replicate_type = None  # Reset replicate_type after assigning sequence
    
    return transcript_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse CDS start and stop positions from a file.")
    parser.add_argument("file_path", type=str, help="Path to the input file with gene sequences")
    parser.add_argument("source_gene_path", type=str, help="Path to the source gene file")
    parser.add_argument("original_gene_path", type=str, help="Path to the original gene file")
    args = parser.parse_args()
    
    RewriteDNA(args.file_path, args.source_gene_path, args.original_gene_path)
