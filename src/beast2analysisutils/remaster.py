import dendropy
from Bio import AlignIO
from xml.etree import ElementTree as ET
from xml.dom import minidom
from datetime import datetime, timedelta
import pandas as pd
import os
import io
import re

def get_state0_newick(nexus_text):
    """
    Given the full text of a NEXUS file, extract the tree (STATE_0 or TREE)
    as a Newick string with integer IDs replaced by leaf names.
    """

    # 1. Extract the TRANSLATE block to build a dictionary:
    #    numeric ID => leaf name (if present)
    translate_block_pattern = re.compile(
        r"(?is)Begin\s+trees\s*;.*?Translate\s*(.*?)\s*;", re.DOTALL
    )
    translate_match = translate_block_pattern.search(nexus_text)
    translation_dict = {}
    if translate_match:
        translate_block = translate_match.group(1)
        # Build a dictionary from the translate block
        for line in translate_block.split("\n"):
            line = line.strip().rstrip(",;")
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                num_id, label = parts
                translation_dict[num_id] = label
    else:
        # No Translate block - will use numeric IDs directly
        pass

    # 2. Extract the tree STATE_0 or TREE = ( ... );
    #    We'll match up to the semicolon
    state0_pattern = re.compile(r"tree\s+STATE_0\s*=\s*(.*?)\s*;", re.DOTALL)
    state0_match = state0_pattern.search(nexus_text)

    if not state0_match:
        # Try TREE instead
        state0_pattern = re.compile(r"tree\s+TREE\s*=\s*(.*?)\s*;", re.DOTALL)
        state0_match = state0_pattern.search(nexus_text)

    if not state0_match:
        raise ValueError(
            "Could not find 'tree STATE_0 = ...;' or 'tree TREE = ...;' in the NEXUS text."
        )

    tree_str = state0_match.group(1).strip()

    # 3. Remove bracketed metadata like [&type="D0",time=0.0]
    tree_str = re.sub(r"\[.*?\]", "", tree_str)

    # 4. Safely replace only integer IDs with leaf names
    #    (leave floating-point numbers alone)
    def replace_labels_safely(s, dictionary):
        result = []
        i = 0
        while i < len(s):
            c = s[i]
            if c.isdigit():
                j = i + 1
                while j < len(s) and s[j].isdigit():
                    j += 1
                numeric_token = s[i:j]
                if j < len(s) and s[j] == ":":
                    if numeric_token in dictionary:
                        result.append(dictionary[numeric_token])
                    else:
                        result.append(numeric_token)
                else:
                    result.append(numeric_token)
                i = j
            else:
                result.append(c)
                i += 1
        return "".join(result)

    tree_str = replace_labels_safely(tree_str, translation_dict)

    # 5. Remove extraneous whitespace, then return
    tree_str = re.sub(r"\s+", "", tree_str)

    # Return a valid Newick with trailing semicolon
    return tree_str + ";"

def extract_remaster_data(alignment_path, tree_path):
    """
    Extracts sequences, leaf dates, and leaf types from ReMASTER output files.

    Args:
        alignment_path (str or None): Path to the Nexus alignment file. If None, sequences will be artificial (empty/placeholder).
        tree_path (str): Path to the Nexus tree file.

    Returns:
        tuple: (sequences, dates, types, newick_tree)
            - sequences (dict): {taxon_id: sequence_string}
            - dates (dict): {taxon_id: date_float_or_string}
            - types (dict): {taxon_id: type_string}
            - newick_tree (str): The tree as a Newick string.
    """
    # 1. Parse Alignment
    sequences = {}
    if alignment_path:
        alignment = AlignIO.read(alignment_path, "nexus")
        for record in alignment:
            sequences[record.id] = str(record.seq).lower()
    else:
        # If no alignment, we will populate dummy sequences based on tree leaves later
        pass

    # 2. Parse Tree and extract metadata
    # We need to handle the potential trailing comma in TRANSLATE block issue mentioned in original script
    # For now, let's try reading directly with dendropy, if it fails, we implement the fix.
    # The original script implemented a fix for 'TRANSLATE' block trailing comma. 
    # Let's incorporate that fix to be safe.
    
    with open(tree_path, "r") as f:
        lines = f.readlines()

    in_translate = False
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("translate"):
            in_translate = True
            continue
        if in_translate:
            if ";" in line:
                prev_idx = i - 1
                while prev_idx > 0 and lines[prev_idx].strip() == "":
                    prev_idx -= 1
                if prev_idx > 0:
                    prev_line = lines[prev_idx].rstrip("\n")
                    if prev_line.rstrip().endswith(","):
                        lines[prev_idx] = prev_line.rstrip(",\n") + "\n"
                in_translate = False
                break
    
    # Use io.StringIO to avoid writing temp file to disk if possible, 
    # but dendropy .get(path=...) or .get(data=...)
    tree_content = "".join(lines)
    tree = dendropy.Tree.get(
        data=tree_content, schema="nexus", preserve_underscores=True
    )
    
    times = {}
    types = {}
    
    for leaf in tree.leaf_node_iter():
        real_label = leaf.taxon.label
        
        # Extract annotations
        # ReMASTER/BEAST2 annotations usually in [&type="...", time=...]
        # dendropy parses these into annotations
        
        if leaf.annotations:
            # Note: The original script accessed 'type' and 'time'.
            # Let's try to get them safely.
            type_val = leaf.annotations.get_value("type")
            time_val = leaf.annotations.get_value("time")
            
            if type_val:
                # Remove curly braces if present (some versions of dendropy/remaster might keep them?)
                # The original script did .replace('{','').replace('}','')
                types[real_label] = str(type_val).replace('{','').replace('}','')
            
            if time_val is not None:
                times[real_label] = float(time_val)
    
    # If sequences dict is empty (no alignment provided), fill with dummy
    if not sequences:
        # Just use "n" as dummy sequence for each taxon found in tree
        # Length doesn't matter much for fixed tree if not used for likelihood, 
        # but BEAST XML usually requires some sequence data.
        # However, for fixed tree, often we don't have sequences? 
        # User requirement: "sequences still need to be placed but they are artifical"
        for taxon in times.keys():
            sequences[taxon] = "n"

    # Also extract the raw Newick tree string for fixed tree templates
    with open(tree_path, "r") as f:
        tree_text = f.read()
    
    # We use the custom parser to get a clean Newick string suitable for TreeParser
    try:
        newick_tree = get_state0_newick(tree_text)
    except Exception as e:
        # If extraction fails (e.g. unexpected format), we might not have a tree string
        # This is okay if we aren't using a fixed tree template, but let's warn or return None
        print(f"Warning: Could not extract Newick tree string: {e}")
        newick_tree = None

    return sequences, times, types, newick_tree

def convert_times_to_dates(times, start_date="2000/01/01"):
    """
    Converts relative simulation times (in years) to date strings based on a start date.

    Args:
        times (dict): {taxon_id: relative_time_float}
        start_date (str): Format "YYYY/MM/DD", defaults to "2000/01/01".

    Returns:
        dict: {taxon_id: date_string} where date_string is "YYYY/MM/DD".
    """
    start_date_dt = datetime.strptime(start_date, "%Y/%m/%d")
    dates = {}
    
    # Logic extracted from user snippet:
    # rel_time is in years, convert to days
    
    for taxon, rel_time in times.items():
        rel_time_days = int(float(rel_time) * 365)
        new_date = start_date_dt + timedelta(days=rel_time_days)
        dates[taxon] = new_date.strftime("%Y/%m/%d")
        
    return dates

def fill_template(template_path, output_path, sequences, dates, types, newick_tree=None):
    """
    Fills the BEAST2 XML template with the provided data.

    Args:
        template_path (str): Path to the template XML file.
        output_path (str): Path to save the generated XML file.
        sequences (dict): {taxon_id: sequence_string}
        dates (dict): {taxon_id: date_value_string} - properly formatted dates e.g. YYYY/MM/DD
        types (dict): {taxon_id: type_value}
        newick_tree (str, optional): The Newick tree string to insert (for fixed tree templates).
    """
    
    # Read template content
    with open(template_path, "r") as f:
        template_content = f.read()

    # Validate date formats
    for taxon, date_str in dates.items():
        try:
            # Check for the primary format we use
            datetime.strptime(date_str, "%Y/%m/%d")
        except ValueError:
            raise ValueError(f"Invalid date format for taxon '{taxon}': '{date_str}'. Expected format YYYY/MM/DD.")

    # 1. Prepare Sequences XML Block
    # Format: <sequence id="seq_{taxon}" spec="Sequence" taxon="{taxon}" totalcount="4" value="{seq}"/>
    seq_lines = []
    # Sort keys for deterministic output
    for taxon in sorted(sequences.keys()):
        seq = sequences[taxon]
        line = f'    <sequence id="seq_{taxon}" spec="Sequence" taxon="{taxon}" totalcount="4" value="{seq}"/>'
        seq_lines.append(line)
    
    sequences_xml_block = "\n".join(seq_lines)
    
    # 2. Prepare Dates String
    # Format: taxon=date,taxon=date...
    date_strings = []
    for taxon in sorted(dates.keys()):
        date_strings.append(f"{taxon}={dates[taxon]}")
    dates_value = ",".join(date_strings)
    
    # 3. Prepare Types String
    # Format: taxon=type,taxon=type...
    type_strings = []
    for taxon in sorted(types.keys()):
        type_strings.append(f"{taxon}={types[taxon]}")
    types_value = ",".join(type_strings)
    
    # 4. Perform Replacements
    # We use simple string replacement as requested for specific placeholders
    # It is safer than parsing XML if the placeholders are just text markers inside tags
    # typically these placeholders are inside value="..." attributes or just as text content.
    
    new_content = template_content
    
    if "INSERTSEQUENCES" in new_content:
        new_content = new_content.replace("INSERTSEQUENCES", sequences_xml_block)
    
    if "INSERTTRAITDATES" in new_content:
        # Only replace if we have headers; otherwise maybe we should remove the block?
        # But commonly we just put the dates.
        new_content = new_content.replace("INSERTTRAITDATES", dates_value)
        
    if "INSERTTRAITTYPES" in new_content:
        new_content = new_content.replace("INSERTTRAITTYPES", types_value)

    if "INSERTNEWICKTREE" in new_content:
        if newick_tree:
            new_content = new_content.replace("INSERTNEWICKTREE", newick_tree)
        else:
            raise ValueError("Template requires INSERTNEWICKTREE but no tree string was provided.")
        
    # Write output
    with open(output_path, "w") as f:
        f.write(new_content)

def generate_xml(template_path, tree_path, output_path, nexus_path=None, fixed_tree=False, start_date="2000/01/01"):
    """
    Wrapper function to generate BEAST2 XML from ReMASTER output.

    Args:
        template_path (str): Path to the XML template.
        tree_path (str): Path to the ReMASTER tree file.
        output_path (str): Path where the final XML will be saved.
        nexus_path (str, optional): Path to the ReMASTER nexus alignment file. Required if fixed_tree is False.
        fixed_tree (bool): Whether to use fixed tree mode.
        start_date (str): Start date for time conversion (YYYY/MM/DD). Defaults to "2000/01/01".
    """
    
    # Validation
    if not fixed_tree and nexus_path is None:
        raise ValueError("nexus_path is required when fixed_tree is False.")
    
    # Extract data
    # extract_remaster_data now handles None nexus_path by generating dummy sequences
    sequences, times, types, newick_tree = extract_remaster_data(nexus_path, tree_path)
    
    # Print Info
    print("--- Remaster Wrapper Info ---")
    if nexus_path:
        print(f"Number of sequences in nexus file: {len(sequences)}")
    else:
        print(f"Number of sequences (artificial): {len(sequences)}")
        
    print(f"Number of leaves in tree: {len(times)}")
    
    # Verification of consistency
    if len(sequences) != len(times):
        print(f"WARNING: Mismatch between number of sequences ({len(sequences)}) and tree leaves ({len(times)}).")
    
    # Iterate type trait stats
    trait_counts = {}
    for taxon, trait in types.items():
        trait_counts[trait] = trait_counts.get(trait, 0) + 1
    
    print("\nLeaves per type trait:")
    for trait, count in trait_counts.items():
        print(f"  {trait}: {count}")
        
    # Date info
    if not fixed_tree:
        if times:
            min_time = min(times.values())
            max_time = max(times.values())
            # Convert to dates for display
            # times are relative years from root usually? Or simply floats.
            # convert_times_to_dates logic: start_date + relative_time
            # Assuming 'times' from extract_remaster_data are appropriate for this conversion.
            # Let's show the converted range.
            
            # Temporary internal conversion for display
            temp_dates = convert_times_to_dates(times, start_date)
            date_values = [datetime.strptime(d, "%Y/%m/%d") for d in temp_dates.values()]
            min_date = min(date_values).strftime("%Y/%m/%d")
            max_date = max(date_values).strftime("%Y/%m/%d")
            
            print(f"\nDate Range: {min_date} to {max_date}")
        else:
            print("\nDate Range: N/A (no times found)")
    else:
        print("\nDate Range: Ignored (Fixed Tree mode)")

    print("-----------------------------")

    # Process Data
    real_dates = convert_times_to_dates(times, start_date)
    
    # Generate XML
    fill_template(
        template_path, 
        output_path, 
        sequences, 
        real_dates, 
        types, 
        newick_tree=newick_tree if fixed_tree else None
    )
    print(f"\nSuccessfully generated XML at: {output_path}")
