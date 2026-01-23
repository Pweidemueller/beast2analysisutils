import os
import pytest
from beast2analysisutils.remaster import extract_remaster_data, fill_template, convert_times_to_dates, get_state0_newick

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
# Templates are now in the package source
PACKAGE_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "beast2analysisutils")
TEMPLATE_PATH = os.path.join(PACKAGE_DIR, "data", "templates", "Template.xml")
FIXED_TEMPLATE_PATH = os.path.join(PACKAGE_DIR, "data", "templates", "Template_fixedtree.xml")

ALIGNMENT_PATH = os.path.join(DATA_DIR, "1_2_simulation.nexus")
TREE_PATH = os.path.join(DATA_DIR, "1_2_simulation.trees")

MOCK_TREE = """
#NEXUS
begin trees;
    Translate
        1 taxonA,
        2 taxonB;
    tree STATE_0 = [&R] (1:1.0,2:1.0):0.0;
end;
"""

def test_get_state0_newick():
    newick = get_state0_newick(MOCK_TREE)
    assert newick.strip() == "(taxonA:1.0,taxonB:1.0):0.0;"

@pytest.mark.skipif(not os.path.exists(ALIGNMENT_PATH), reason="Sample data not found")
def test_remaster_workflow(tmp_path):
    # 1. Extract data (including tree now)
    sequences, times, types, newick_tree = extract_remaster_data(ALIGNMENT_PATH, TREE_PATH)
    
    assert len(sequences) > 0
    assert len(times) > 0
    assert len(types) > 0
    # newick_tree might be None or a string depending on implementation details of existing function if tree file is not perfect
    # but with sample data it should be present.
    # Actually, in the updated code, extract_remaster_data returns 4 values.
    
    # Check consistency
    assert len(sequences) == len(times)
    
    # Convert times to dates
    dates = convert_times_to_dates(times, start_date="2000/01/01")
    assert len(dates) == len(times)
    
    # 2. Fill template (Standard)
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

@pytest.mark.skipif(not os.path.exists(ALIGNMENT_PATH), reason="Sample data not found")
def test_remaster_fixedtree_workflow(tmp_path):
    sequences, times, types, newick_tree = extract_remaster_data(ALIGNMENT_PATH, TREE_PATH)
    dates = convert_times_to_dates(times, start_date="2000/01/01")
    
    assert newick_tree is not None
    assert newick_tree.strip().endswith(";")
    
    output_xml = tmp_path / "output_fixed.xml"
    
    # If Template_fixedtree.xml exists, use it, otherwise mock it for logic test
    template_path = FIXED_TEMPLATE_PATH
    if not os.path.exists(FIXED_TEMPLATE_PATH):
        # Create a mock fixed template
        mock_template = tmp_path / "mock_fixed_template.xml"
        with open(mock_template, "w") as f:
            f.write('<beast><data>INSERTSEQUENCES</data><tree>INSERTNEWICKTREE</tree><trait>INSERTTRAITTYPES</trait></beast>')
        template_path = str(mock_template)

    fill_template(template_path, str(output_xml), sequences, dates, types, newick_tree=newick_tree)
    
    assert output_xml.exists()
    with open(output_xml, "r") as f:
        content = f.read()
        
    assert "INSERTNEWICKTREE" not in content
    # Check tree string is inserted
    # We strip whitespace from expected tree as fill_template might not modify whitespace but XML might be messy?
    # Actually fill_template does replace directly.
    assert newick_tree in content

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

def test_generate_xml_wrapper_standard(tmp_path):
    from beast2analysisutils.remaster import generate_xml
    
    output_xml = tmp_path / "wrapper_output.xml"
    if not os.path.exists(ALIGNMENT_PATH):
        pytest.skip("Alignment file needed for standard test")
        
    generate_xml(
        template_path=TEMPLATE_PATH,
        tree_path=TREE_PATH,
        output_path=str(output_xml),
        alignment_path=ALIGNMENT_PATH,
        fixed_tree=False
    )
    
    assert output_xml.exists()
    # Check basics
    with open(output_xml, "r") as f:
        content = f.read()
    assert "INSERTSEQUENCES" not in content

def test_generate_xml_wrapper_fixed_tree(tmp_path):
    from beast2analysisutils.remaster import generate_xml
    
    output_xml = tmp_path / "wrapper_output_fixed.xml"
    
    # Mock fixed template if not exists
    template_path = FIXED_TEMPLATE_PATH
    if not os.path.exists(FIXED_TEMPLATE_PATH):
        mock_template = tmp_path / "mock_fixed_template.xml"
        with open(mock_template, "w") as f:
            f.write('<beast><data>INSERTSEQUENCES</data><tree>INSERTNEWICKTREE</tree><trait>INSERTTRAITTYPES</trait></beast>')
        template_path = str(mock_template)

    # Run without alignment (optional for fixed tree)
    generate_xml(
        template_path=template_path,
        tree_path=TREE_PATH,
        output_path=str(output_xml),
        alignment_path=None,
        fixed_tree=True
    )
    
    assert output_xml.exists()
    with open(output_xml, "r") as f:
        content = f.read()
    # Should have artificial sequences "n"
    assert 'value="n"' in content or 'value="N"' in content

def test_generate_xml_validation():
    from beast2analysisutils.remaster import generate_xml
    
    # Should raise error if fixed_tree=False but no alignment_path
    with pytest.raises(ValueError, match="alignment_path is required"):
        generate_xml("tpl", "tree", "out", alignment_path=None, fixed_tree=False)
