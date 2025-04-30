"""
# utils.py

Miscellaneous utility functions used by this project
"""


def regularize_name(name: str):
    """Regularizes a name for comparisons, making it lowercase and stripping spaces

    Parameters
    ----------
    name : str
        The name, e.g. "Open Babel"

    Returns
    -------
    str
        The regularized name, e.g. "openbabel"
    """
    return name.lower().replace(" ", "")
