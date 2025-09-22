import json
import os

def parse_fasta_universal(file_path, make_uppercase=False):
    """
    FASTA parser that handles all provided formats by separating the header from the sequence.
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            # Split the entire file into records based on the '>' character
            records = f.read().strip().split('>')

            for record_block in records:
                # Skip any empty records that might result from splitting (e.g., the first one)
                if not record_block.strip():
                    continue

                # Find the position of the first newline character
                try:
                    first_newline_pos = record_block.index('\n')
                    # The part of the string *after* the first newline is the sequence
                    sequence_part = record_block[first_newline_pos + 1:]
                except ValueError:
                    # This handles cases where a record might be just a header with no sequence
                    sequence_part = ""

                # Remove all whitespace (like newlines) from the sequence part to join it
                cleaned_sequence = "".join(sequence_part.split())

                # Special handling for mouse data to make it uppercase
                if make_uppercase:
                    cleaned_sequence = cleaned_sequence.upper()

                # Only add the sequence if it's not empty
                if cleaned_sequence:
                    sequences.append(cleaned_sequence)

    except FileNotFoundError:
        print(f"Warning: File not found at {file_path}. Skipping.")
        return None # Return None to indicate failure

    return sequences

def process_react_numeric(file_path):
    """Reads a file with one numeric reactivity value per line, skipping the header."""
    react_values = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        if lines and "reacts" in lines[0].lower():
            lines = lines[1:] # Skip header
        for line in lines:
            if line.strip():
                react_values.append(int(line.strip()))
    return react_values

def process_react_binary(file_path):
    """Reads a file with binary 'Y'/'N' reactivity and converts to 1/0."""
    react_values = []
    with open(file_path, 'r') as f:
        for line in f:
            val = line.strip().upper()
            if val == 'Y':
                # This will be set to 4 (not 1) so it doesn't get affected by downstream data filtering:
                # 0:   non-specific
                # 1-3: mildly specific (removed from dataset)
                # >3:  specific
                react_values.append(4)
            elif val == 'N':
                react_values.append(0)
    return react_values

# --- Main Execution ---

all_data = []
data_directory = '.'
prefixes = ['flu', 'gut_hiv', 'nat_cntrl', 'nat_hiv', 'plos_hiv', 'mouse']

for prefix in prefixes:
    print(f"Processing {prefix} data...")

    # Determine file extensions for sequence files
    h_ext = '.dat' if prefix == 'mouse' else '.txt'
    l_ext = '.dat' if prefix == 'mouse' else '.txt'
    
    # Uppercase flag for mouse data
    should_uppercase = (prefix == 'mouse')

    # Define file paths
    heavy_chain_file = os.path.join(data_directory, f'{prefix}_fastaH{h_ext}')
    light_chain_file = os.path.join(data_directory, f'{prefix}_fastaL{l_ext}')

    # Parse the heavy and light chain sequences
    heavy_chains = parse_fasta_universal(heavy_chain_file, make_uppercase=should_uppercase)
    light_chains = parse_fasta_universal(light_chain_file, make_uppercase=should_uppercase)

    # If parsing failed for any file, skip this entire prefix
    if heavy_chains is None or light_chains is None:
        continue

    # Determine which reactivity parser to use
    if prefix in ['mouse', 'plos_hiv']:
        reactivity_file = os.path.join(data_directory, f'{prefix}_YN.txt')
        reactivity = process_react_binary(reactivity_file)
    else:
        reactivity_file = os.path.join(data_directory, f'{prefix}_NumReact.txt')
        reactivity = process_react_numeric(reactivity_file)

    # Check for data consistency
    if len(heavy_chains) != len(light_chains) or len(heavy_chains) != len(reactivity):
        print(f"Warning: Mismatch in record counts for prefix '{prefix}'. Skipping.")
        print(f"  Heavy chains: {len(heavy_chains)}, Light chains: {len(light_chains)}, Reactivity values: {len(reactivity)}")
        continue

    # Combine the data for the current prefix
    for h_chain, l_chain, react_val in zip(heavy_chains, light_chains, reactivity):
        all_data.append({
            "type": prefix,
            "NumReact": react_val,
            "DNA_fastaH_raw": h_chain,
            "DNA_fastaL_raw": l_chain
        })

# --- Save the consolidated data ---
output_filename = 'consolidated_data_boughter_2020.json'
with open(output_filename, 'w') as f:
    json.dump(all_data, f, indent=2)

print(f"\nSuccessfully consolidated all data into '{output_filename}'")
print(f"Total records processed: {len(all_data)}")