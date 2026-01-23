import dendropy
from Bio import AlignIO
from xml.etree import ElementTree as ET
from xml.dom import minidom
from datetime import datetime, timedelta
import pandas as pd
import os
import io

def extract_remaster_data(alignment_path, tree_path):
    """
    Extracts sequences, leaf dates, and leaf types from ReMASTER output files.

    Args:
        alignment_path (str): Path to the Nexus alignment file.
        tree_path (str): Path to the Nexus tree file.

    Returns:
        tuple: (sequences, dates, types)
            - sequences (dict): {taxon_id: sequence_string}
            - dates (dict): {taxon_id: date_float_or_string}
            - types (dict): {taxon_id: type_string}
    """
    # 1. Parse Alignment
    alignment = AlignIO.read(alignment_path, "nexus")
    sequences = {}
    for record in alignment:
        sequences[record.id] = str(record.seq).lower()

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
        
    return sequences, times, types

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

def fill_template(template_path, output_path, sequences, dates, types):
    """
    Fills the BEAST2 XML template with the provided data.

    Args:
        template_path (str): Path to the template XML file.
        output_path (str): Path to save the generated XML file.
        sequences (dict): {taxon_id: sequence_string}
        dates (dict): {taxon_id: date_value_string} - properly formatted dates e.g. YYYY/MM/DD
        types (dict): {taxon_id: type_value}
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
        new_content = new_content.replace("INSERTTRAITDATES", dates_value)
        
    if "INSERTTRAITTYPES" in new_content:
        new_content = new_content.replace("INSERTTRAITTYPES", types_value)
        
    # Write output
    with open(output_path, "w") as f:
        f.write(new_content)
