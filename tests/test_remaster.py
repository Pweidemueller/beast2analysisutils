import os
import pytest
from beast2analysisutils.remaster import extract_remaster_data, fill_template, convert_times_to_dates

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "temporary_files")
TEMPLATE_PATH = os.path.join(DATA_DIR, "Template.xml")
ALIGNMENT_PATH = os.path.join(DATA_DIR, "1_2_simulation.nexus")
TREE_PATH = os.path.join(DATA_DIR, "1_2_simulation.trees")

@pytest.mark.skipif(not os.path.exists(ALIGNMENT_PATH), reason="Sample data not found")
def test_remaster_workflow(tmp_path):
    # 1. Extract data
    sequences, times, types = extract_remaster_data(ALIGNMENT_PATH, TREE_PATH)
    
    assert len(sequences) > 0
    assert len(times) > 0
    assert len(types) > 0
    
    # Check consistency
    assert len(sequences) == len(times)
    
    # Convert times to dates
    dates = convert_times_to_dates(times, start_date="2000/01/01")
    assert len(dates) == len(times)
    
    # 2. Fill template
    output_xml = tmp_path / "output.xml"
    fill_template(TEMPLATE_PATH, str(output_xml), sequences, dates, types)
    
    assert output_xml.exists()
    
    with open(output_xml, "r") as f:
        content = f.read()
        
    # Check placeholders are gone
    assert "INSERTSEQUENCES" not in content
    assert "INSERTTRAITDATES" not in content
    assert "INSERTTRAITTYPES" not in content
    
    # Check some content was inserted
    # Pick a random taxon from sequences
    taxon = next(iter(sequences))
    assert f'taxon="{taxon}"' in content
    assert f'{taxon}=' in content # in dates or types string

def test_fill_template_invalid_date():
    """Test that fill_template raises ValueError for invalid date formats."""
    # Create dummy data
    sequences = {"t1": "ACGT"}
    # Invalid date format
    dates = {"t1": "2020-01-01"} # Expects YYYY/MM/DD
    types = {"t1": "type1"}
    
    # Use a dummy template file
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp:
        tmp.write("template content")
        tmp_path = tmp.name
    
    try:
        with pytest.raises(ValueError, match="Invalid date format"):
            fill_template(tmp_path, "output.xml", sequences, dates, types)
    finally:
        os.remove(tmp_path)
