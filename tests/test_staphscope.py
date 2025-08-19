import pytest
from pathlib import Path
from staphscope.core import check_environment

def test_environment_check():
    """Test that environment check works"""
    tools = check_environment()
    assert isinstance(tools, dict)
    assert 'python3' in tools

# Add more tests as needed
