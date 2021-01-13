import os
from pathlib import Path
from typing import Sequence, Union


def get_test_data(folder: str, patterns: Union[str, Sequence[str]]):
    testdata_path = os.getenv("RAVENPY_TESTDATA_PATH")
    if not testdata_path:
        raise RuntimeError("RAVENPY_TESTDATA_PATH env variable is not set")
    testdata_path = Path(testdata_path)
    patterns = [patterns] if isinstance(patterns, str) else patterns
    return [
        p
        for pat in patterns
        for p in (testdata_path / folder).glob(pat)
        if p.suffix != ".md5"
    ]
