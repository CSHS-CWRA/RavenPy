"""
Tools for searching for and acquiring test data.
"""
import logging
import os
import re
from hashlib import md5
from pathlib import Path
from typing import List, Optional, Sequence, Union
from urllib.error import HTTPError
from urllib.parse import urljoin
from urllib.request import urlretrieve

import requests
from xarray import Dataset
from xarray import open_dataset as _open_dataset

_default_cache_dir = Path.home() / ".raven_testing_data"

LOGGER = logging.getLogger("RAVEN")

__all__ = ["get_local_testdata", "get_file", "open_dataset", "query_folder"]


def file_md5_checksum(fname):
    hash_md5 = md5()
    with open(fname, "rb") as f:
        hash_md5.update(f.read())
    return hash_md5.hexdigest()


def get_local_testdata(pattern: str) -> Union[Path, List[Path]]:
    """Gather testdata from a local folder.

    Return files matching `pattern` in the local test data repo
    located at `RAVENPY_TESTDATA_PATH` (which must be set).

    Parameters
    ----------
    pattern: str
      Glob pattern, which must include the folder.

    Returns
    -------
    Union[Path, List[Path]]
    """
    testdata_path = os.getenv("RAVENPY_TESTDATA_PATH")
    if not testdata_path:
        raise RuntimeError("RAVENPY_TESTDATA_PATH env variable is not set")
    testdata_path = Path(testdata_path)
    if not testdata_path.exists():
        raise RuntimeError(f"{testdata_path} does not exists")
    paths = [path for path in testdata_path.glob(pattern) if path.suffix != ".md5"]
    if not paths:
        raise RuntimeError(f"No data found for {pattern} at {testdata_path}")
    # Return item directly when singleton, for convenience
    return paths[0] if len(paths) == 1 else paths


def _get(
    fullname: Path,
    github_url: str,
    branch: str,
    suffix: str,
    cache_dir: Path,
) -> Path:
    cache_dir = cache_dir.absolute()
    local_file = cache_dir / branch / fullname
    md5name = fullname.with_suffix("{}.md5".format(suffix))
    md5file = cache_dir / branch / md5name

    if not local_file.is_file():
        # This will always leave this directory on disk.
        # We may want to add an option to remove it.
        local_file.parent.mkdir(parents=True, exist_ok=True)

        url = "/".join((github_url, "raw", branch, fullname.as_posix()))
        LOGGER.info("Fetching remote file: %s" % fullname.as_posix())
        urlretrieve(url, local_file)
        try:
            url = "/".join((github_url, "raw", branch, md5name.as_posix()))
            LOGGER.info("Fetching remote file md5: %s" % md5name.as_posix())
            urlretrieve(url, md5file)
        except HTTPError as e:
            msg = f"{md5name.as_posix()} not found. Aborting file retrieval."
            local_file.unlink()
            raise FileNotFoundError(msg) from e

        localmd5 = file_md5_checksum(local_file)
        try:
            with open(md5file) as f:
                remotemd5 = f.read()
            if localmd5 != remotemd5:
                local_file.unlink()
                msg = """
                    MD5 checksum does not match, try downloading dataset again.
                    """
                raise OSError(msg)
        except OSError as e:
            LOGGER.error(e)
    return local_file


# idea copied from xclim that borrowed it from xarray that was borrowed from Seaborn
def get_file(
    name: Union[str, Sequence[str]],
    github_url: str = "https://github.com/Ouranosinc/raven-testdata",
    branch: str = "master",
    cache_dir: Path = _default_cache_dir,
) -> Union[Path, List[Path]]:
    """
    Return a file from an online GitHub-like repository.
    If a local copy is found then always use that to avoid network traffic.

    Parameters
    ----------
    name : Union[str, Sequence[str]]
        Name of the file or list/tuple of names of files containing the dataset(s) including suffixes.
    github_url : str
        URL to Github repository where the data is stored.
    branch : str, optional
        For GitHub-hosted files, the branch to download from.
    cache_dir : Path
        The directory in which to search for and write cached data.

    Returns
    -------
    Union[Path, List[Path]]
    """
    if isinstance(name, str):
        name = [name]

    files = list()
    for n in name:
        fullname = Path(n)
        suffix = fullname.suffix
        files.append(
            _get(
                fullname=fullname,
                github_url=github_url,
                branch=branch,
                suffix=suffix,
                cache_dir=cache_dir,
            )
        )
    if len(files) == 1:
        return files[0]
    return files


# Credits to Anselme  https://stackoverflow.com/a/62003257/7322852 (CC-BY-SA 4.0)
def query_folder(
    folder: Optional[str] = None,
    pattern: Optional[str] = None,
    github_url: str = "https://github.com/Ouranosinc/raven-testdata",
    branch: str = "master",
) -> List[str]:
    """
    Lists the files available for retrieval from a remote git repository with get_file.
    If provided a folder name, will perform a globbing-like filtering operation for parent folders.

    Parameters
    ----------
    folder : str, optional
        Relative pathname of the sub-folder from the top-level.
    pattern : str, optional
        Regex pattern to identify a file.
    github_url : str
        URL to Github repository where the data is stored.
    branch : str, optional
        For GitHub-hosted files, the branch to download from.

    Returns
    -------
    List[str]
    """
    repo_name = github_url.strip("https://github.com/")

    url = f"https://api.github.com/repos/{repo_name}/git/trees/{branch}?recursive=1"
    r = requests.get(url)
    res = r.json()

    md5_files = [f["path"] for f in res["tree"] if f["path"].endswith(".md5")]
    if folder:
        folder = "/".join("/".split(folder)) if "/" in folder else folder
        md5_files = [f for f in md5_files if folder in Path(f).parent.as_posix()]
    files = [re.sub(".md5$", "", f) for f in md5_files]

    if pattern:
        regex = re.compile(pattern)
        files = [string for string in files if re.search(regex, string)]

    return files


# idea copied from xclim that borrowed it from xarray that was borrowed from Seaborn
def open_dataset(
    name: str,
    suffix: Optional[str] = None,
    dap_url: Optional[str] = None,
    github_url: str = "https://github.com/Ouranosinc/raven-testdata",
    branch: str = "master",
    cache: bool = True,
    cache_dir: Path = _default_cache_dir,
    **kwds,
) -> Dataset:
    """Open a dataset from the online GitHub-like repository.

    If a local copy is found then always use that to avoid network traffic.

    Parameters
    ----------
    name: str
      Name of the file containing the dataset. If no suffix is given, assumed to be netCDF ('.nc' is appended).
    suffix: str, optional
      If no suffix is given, assumed to be netCDF ('.nc' is appended). For no suffix, set "".
    dap_url: str, optional
      URL to OPeNDAP folder where the data is stored. If supplied, supersedes github_url.
    github_url: str
      URL to Github repository where the data is stored.
    branch: str, optional
      For GitHub-hosted files, the branch to download from.
    cache_dir: Path
      The directory in which to search for and write cached data.
    cache: bool
      If True, then cache data locally for use on subsequent calls.
    kwds: dict, optional
      For NetCDF files, keywords passed to xarray.open_dataset.

    Returns
    -------
    Union[Dataset, Path]

    See Also
    --------
    xarray.open_dataset
    """
    name = Path(name)
    if suffix is None:
        suffix = ".nc"
    fullname = name.with_suffix(suffix)

    if dap_url is not None:
        dap_file = urljoin(dap_url, str(name))
        try:
            ds = _open_dataset(dap_file, **kwds)
            return ds
        except OSError:
            msg = "OPeNDAP file not read. Verify that service is available."
            LOGGER.error(msg)
            raise

    local_file = _get(
        fullname=fullname,
        github_url=github_url,
        branch=branch,
        suffix=suffix,
        cache_dir=cache_dir,
    )

    try:
        ds = _open_dataset(local_file, **kwds)
        if not cache:
            ds = ds.load()
            local_file.unlink()
        return ds
    except OSError:
        raise
