import os

def _convert_to_json_serializable(obj):
    """Convert Output objects to nested dictionaries of floats for JSON serialization."""
    if hasattr(obj, '__dict__'):
        # Handle Output objects by converting to dict
        result = {}
        for key, value in dict(obj).items():
            result[key] = _convert_to_json_serializable(value)
        return result
    elif isinstance(obj, (int, float)):
        return float(obj)  # Ensure all numbers are floats
    elif isinstance(obj, str):
        return obj
    elif obj is None:
        return None
    elif isinstance(obj, (list, tuple)):
        return [_convert_to_json_serializable(item) for item in obj]
    elif isinstance(obj, dict):
        return {k: _convert_to_json_serializable(v) for k, v in obj.items()}
    else:
        # For other types, try to convert to string
        return str(obj)

def get_filepath(filename: str) -> str:
    """
    Standardize filepath for test files by placing them in a test_outputs directory.
    
    Args:
        filename: Base filename for the test
        
    Returns:
        Full filepath with standardized directory
    """
    test_dir = "test_reference_interstitial"
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    return os.path.join(test_dir, filename)