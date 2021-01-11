from typing import Sequence, Union

from ravenpy import VENV_PATH

TESTDATA_BASE_FOLDER_IN_VENV = "raven-testdata-master"


def get_test_data(folder: str, patterns: Union[str, Sequence[str]]):
    testdata_path = VENV_PATH / TESTDATA_BASE_FOLDER_IN_VENV
    if not testdata_path.exists():
        raise RuntimeError(
            f"""Cannot find testdata in the expected location: {testdata_path}
Please make sure that RavenPy was pip installed with the `--install-option="--with-testdata"` option, or see README for more information.
"""
        )
    patterns = [patterns] if isinstance(patterns, str) else patterns
    return [
        p
        for pat in patterns
        for p in (testdata_path / folder).glob(pat)
        if p.suffix != ".md5"
    ]
