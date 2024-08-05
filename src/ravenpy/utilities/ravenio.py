"""Tools for reading outputs and writing inputs for the Raven executable."""

import re
from collections import OrderedDict
from typing import Any


# TODO: Implement section parser
def parse_configuration(fn) -> dict[str, Any]:
    """
    Parse Raven configuration file.

    Parameters
    ----------
    fn : str or Path
        Path to the configuration file.

    Returns
    -------
    dict
        A dictionary keyed by parameter name.
    """
    main_param = re.compile(r"^:(\w+)\s+([^#]*)")
    # sub_param = re.compile(r"^  :(\w+)\s+([^#]*)")
    out = OrderedDict()
    # cat = None
    with Path(str(fn)).open() as f:
        for line in f.readlines():
            match = main_param.search(line)
            if not match:
                continue

            key, value = match.groups()
            if value:
                values = value.split()
                out[key] = values[0] if len(values) == 1 else values
            else:
                if "List" in key:
                    pass
                elif "Classes" in key:
                    pass
                elif "Profiles" in key:
                    pass
                else:
                    out[key] = True

    return out
