import os
import os.path
import tempfile
import urllib.request
import zipfile
from datetime import datetime
from importlib.metadata import version
from os import PathLike
from pathlib import Path
from textwrap import dedent, indent
from typing import cast

import mudata
import pandas as pd
import pooch
import scanpy as sc
from anndata import AnnData
from mudata import MuData
from pandas.api.types import is_string_dtype
from scanpy import logging

from scirpy.io._convert_anndata import from_airr_cells
from scirpy.io._datastructures import AirrCell
from scirpy.io._io import _infer_locus_from_gene_names, _IOLogger
from scirpy.pp import index_chains
from scirpy.util import _doc_params, _read_to_str, tqdm

HERE = Path(__file__).parent

_FIGSHARE = pooch.create(
    path=pooch.os_cache("scirpy"),
    base_url="doi:10.6084/m9.figshare.22249894.v1",
    version=version("scirpy"),
    version_dev="main",
    env="SCIRPY_DATA_DIR",
    registry={
        "wu2020.h5mu": "md5:ed30d9c1c44cae544f4c080a2451118b",
        "wu2020_3k.h5mu": "md5:12c57c790f8a403751304c9de5a18cbf",
        "maynard2020.h5mu": "md5:da64ac62e3e92c80eaf0e8eef6537ac7",
    },
)
_POOCH_INFO = dedent(
    """\
    .. note::
        Scirpy example datasets are managed through `Pooch <https://github.com/fatiando/pooch>`_.

        By default, the dataset will be downloaded into your operating system's default
        cache directory (See :func:`pooch.os_cache` for more details). If it has already been
        downloaded, it will be retrieved from the cache.

        You can override the default cache dir by setting the `SCIRPY_DATA_DIR` environment variable
        to a path of your preference.
    """
)


@_doc_params(
    processing_code=indent(_read_to_str(HERE / "_processing_scripts/wu2020.py"), " " * 8),
    pooch_info=_POOCH_INFO,
)
def wu2020() -> MuData:
    """\
    Return the dataset from :cite:`Wu2020` as MuData object.

    140k cells, of which 100k have TCRs.

    {pooch_info}

    This is how the dataset was processed:

    .. code-block:: python

        {processing_code}
    """
    fname = cast(PathLike, _FIGSHARE.fetch("wu2020.h5mu", progressbar=True))
    return mudata.read_h5mu(fname)


@_doc_params(
    processing_code=indent(_read_to_str(HERE / "_processing_scripts/wu2020_3k.py"), " " * 8),
    pooch_info=_POOCH_INFO,
)
def wu2020_3k() -> MuData:
    """\
    Return the dataset from :cite:`Wu2020` as AnnData object, downsampled
    to 3000 TCR-containing cells.

    {pooch_info}

    This is how the dataset was processed:

    .. code-block:: python

        {processing_code}
    """
    fname = cast(PathLike, _FIGSHARE.fetch("wu2020_3k.h5mu", progressbar=True))
    return mudata.read_h5mu(fname)


@_doc_params(
    processing_code=indent(_read_to_str(HERE / "_processing_scripts/maynard2020.py"), " " * 8),
    pooch_info=_POOCH_INFO,
)
def maynard2020() -> MuData:
    """\
    Return the dataset from :cite:`Maynard2020` as AnnData object.

    21k cells from NSCLC profiled with Smart-seq2, of which 3,500 have :term:`TCRs<TCR>`
    and 1,500 have :term:`BCRs<BCR>`.

    {pooch_info}

    The raw FASTQ files have been obtained from `PRJNA591860 <https://www.ebi.ac.uk/ena/browser/view/PRJNA591860>`__
    and processed using the nf-core `RNA-seq pipeline <https://github.com/nf-core/rnaseq>`_ to obtain
    gene expression and TraCeR/BraCeR to reconstruct receptors.

    The processed files have been imported and transformed into an :class:`anndata.AnnData`
    object using the following script:

    .. code-block:: python

        {processing_code}
    """
    fname = cast(PathLike, _FIGSHARE.fetch("maynard2020.h5mu", progressbar=True))
    return mudata.read_h5mu(fname)


def vdjdb(cached: bool = True, *, cache_path="data/vdjdb.h5ad") -> AnnData:
    """\
    Download VDJdb and process it into an AnnData object.

    `VDJdb <https://vdjdb.cdr3.net/>`_ :cite:`vdjdb` is a curated database of
    T-cell receptor (TCR) sequences with known antigen specificities.

    Parameters
    ----------
    cached
        If `True`, attempt to read from the `data` directory before downloading
    cache_path
        Location where the h5ad object will be saved

    Returns
    -------
    An anndata object containing all entries from VDJDB in `obsm["airr"]`.
    Each entry is represented as if it was a cell, but without gene expression.
    Metadata is stored in `adata.uns["DB"]`.
    """
    if cached:
        try:
            return sc.read_h5ad(cache_path)
        except OSError:
            pass

    logging.info("Downloading latest version of VDJDB")
    with urllib.request.urlopen(
        "https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt"
    ) as url:
        latest_versions = url.read().decode().split()
    url = latest_versions[0]

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        urllib.request.urlretrieve(url, d / "vdjdb.tar.gz")
        with zipfile.ZipFile(d / "vdjdb.tar.gz") as zf:
            zf.extractall(d)
        df = pd.read_csv(d / "vdjdb_full.txt", sep="\t", low_memory=False)

    tcr_cells = []
    for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing VDJDB entries"):
        cell = AirrCell(cell_id=str(idx))
        if not pd.isnull(row["cdr3.alpha"]):
            alpha_chain = AirrCell.empty_chain_dict()
            alpha_chain.update(
                {
                    "locus": "TRA",
                    "junction_aa": row["cdr3.alpha"],
                    "v_call": row["v.alpha"],
                    "j_call": row["j.alpha"],
                    "consensus_count": 0,
                    "productive": True,
                }
            )
            cell.add_chain(alpha_chain)

        if not pd.isnull(row["cdr3.beta"]):
            beta_chain = AirrCell.empty_chain_dict()
            beta_chain.update(
                {
                    "locus": "TRB",
                    "junction_aa": row["cdr3.beta"],
                    "v_call": row["v.beta"],
                    "d_call": row["d.beta"],
                    "j_call": row["j.beta"],
                    "consensus_count": 0,
                    "productive": True,
                }
            )
            cell.add_chain(beta_chain)

        INCLUDE_CELL_METADATA_FIELDS = [
            "species",
            "mhc.a",
            "mhc.b",
            "mhc.class",
            "antigen.epitope",
            "antigen.gene",
            "antigen.species",
            "reference.id",
            "method.identification",
            "method.frequency",
            "method.singlecell",
            "method.sequencing",
            "method.verification",
            "meta.study.id",
            "meta.cell.subset",
            "meta.subject.cohort",
            "meta.subject.id",
            "meta.replica.id",
            "meta.clone.id",
            "meta.epitope.id",
            "meta.tissue",
            "meta.donor.MHC",
            "meta.donor.MHC.method",
            "meta.structure.id",
            "vdjdb.score",
        ]
        for f in INCLUDE_CELL_METADATA_FIELDS:
            cell[f] = row[f]
        tcr_cells.append(cell)

    logging.info("Converting to AnnData object")
    adata = from_airr_cells(tcr_cells)
    index_chains(adata)

    adata.uns["DB"] = {"name": "VDJDB", "date_downloaded": datetime.now().isoformat()}

    # store cache
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    adata.write_h5ad(cast(os.PathLike, cache_path))

    return adata


def iedb(cached: bool = True, *, cache_path="data/iedb.h5ad") -> AnnData:
    """\
    Download IEBD v3 and process it into an AnnData object.

    :cite:`iedb` is a curated database of
    T-cell receptor (TCR) sequences with known antigen specificities.

    Parameters
    ----------
    cached
        If `True`, attempt to read from the `data` directory before downloading
    cache_path
        Location where the h5ad object will be saved

    Returns
    -------
    An anndata object containing all entries from IEDB in `obsm["airr"]`.
    Each entry is represented as if it was a cell, but without gene expression.
    Metadata is stored in `adata.uns["DB"]`.
    """
    if cached:
        try:
            return sc.read_h5ad(cache_path)
        except OSError:
            pass
    logger = _IOLogger()

    url = "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        urllib.request.urlretrieve(url, d / "receptor_full_v3.zip")
        with zipfile.ZipFile(d / "receptor_full_v3.zip") as zf:
            zf.extractall(d)
        iedb_df = pd.concat(
            [
                pd.read_csv(
                    d / x,
                    index_col=None,
                    header=[0, 1],
                    sep=",",
                    na_values=["None"],
                    true_values=["True"],
                    low_memory=False,
                )
                for x in ["tcr_full_v3.csv", "bcr_full_v3.csv"]
            ]
        ).reset_index(drop=True)

    # deal with multiindex (join by " " to single index). This is how columns were
    # named before an IEDB update in 2023-04
    iedb_df.columns = iedb_df.columns.to_series().apply(lambda x: " ".join(x))
    iedb_df = iedb_df.drop_duplicates()
    iedb_df = iedb_df.drop_duplicates(
        [
            "Chain 1 CDR3 Curated",
            "Chain 2 CDR3 Curated",
            "Epitope Source Organism",
            "Epitope Source Molecule",
            "Receptor Reference Name",
            "Chain 1 Type",
            "Chain 2 Type",
        ]
    )

    def replace_curated(input_1, input_2):
        """If no curated CDR3 sequence or V/D/J gene is available, use the calculated one"""
        for input in [input_1, input_2]:
            calculated = input.replace("Curated", "Calculated")
            iedb_df.loc[iedb_df[input].isna(), input] = iedb_df[calculated]
            if is_string_dtype(iedb_df[input]):
                iedb_df[input] = iedb_df[input].str.upper()

    def _get_junction_aa(row, chain="Chain 1"):
        """Compute the "junction_aa" sequence from the protein sequence based on the Start/End indices.

        In the database, only CDR3 (without initial C and terminal W/F) is available.
        """
        try:
            # if the start end is missing (NaN)
            start = int(row[f"{chain} CDR3 Start Curated"])
            end = int(row[f"{chain} CDR3 End Curated"])
        except ValueError:
            return pd.NA
        protein = row[f"{chain} Protein Sequence"]
        if pd.isnull(protein):
            # if a start/end is specified (possibly a "curated" CDR3 is available), but no reference protein sequence given
            return pd.NA

        # +1 at each end for "junction" vs. "cdr3", start one earlier because of different indexing
        junction = protein[start - 2 : end + 1]

        # It turns out that this does not always match - but in the inspected cases, the "curated" CDR3 was
        # incorrectly including the C/F|W.
        # cdr3 = junction[1:-1]
        # print(junction, cdr3, row["Chain 1 CDR3 Curated"])

        return junction

    replace_curated("Chain 1 CDR3 Curated", "Chain 2 CDR3 Curated")
    replace_curated("Chain 1 CDR3 Start Curated", "Chain 2 CDR3 Start Curated")
    replace_curated("Chain 1 CDR3 End Curated", "Chain 2 CDR3 End Curated")
    replace_curated("Chain 1 Curated V Gene", "Chain 2 Curated V Gene")
    replace_curated("Chain 1 Curated D Gene", "Chain 2 Curated D Gene")
    replace_curated("Chain 1 Curated J Gene", "Chain 2 Curated J Gene")

    iedb_df["Chain 1 junction_aa"] = iedb_df.apply(_get_junction_aa, axis=1, chain="Chain 1")
    iedb_df["Chain 2 junction_aa"] = iedb_df.apply(_get_junction_aa, axis=1, chain="Chain 2")
    iedb_df["cell_id"] = iedb_df.reset_index(drop=True).index

    accepted_chains = ["alpha", "beta", "heavy", "light", "gamma", "delta"]
    iedb_df = iedb_df[(iedb_df["Chain 1 Type"].isin(accepted_chains)) & (iedb_df["Chain 2 Type"].isin(accepted_chains))]

    receptor_dict = {
        "alpha": "TRA",
        "beta": "TRB",
        "heavy": "IGH",
        # IEDB does not distinguish between lambda and kappa
        "light": None,
        "gamma": "TRG",
        "delta": "TRD",
    }

    tcr_cells = []
    for _, row in iedb_df.iterrows():
        cell = AirrCell(cell_id=row["cell_id"], logger=logger)
        chain1 = AirrCell.empty_chain_dict()
        chain2 = AirrCell.empty_chain_dict()
        cell["Receptor IEDB Receptor ID"] = row["Receptor IEDB Receptor ID"]
        cell["Epitope Source Molecule"] = row["Epitope Source Molecule"]
        cell["Epitope Source Organism"] = row["Epitope Source Organism"]
        cell["Receptor Reference Namee"] = row["Receptor Reference Name"]
        cell["Reference IEDB IRI"] = row["Reference IEDB IRI"]
        cell["Epitope IEDB IRI"] = row["Epitope IEDB IRI"]
        chain1.update(
            {
                "locus": receptor_dict[row["Chain 1 Type"]],
                "cdr3_aa": row["Chain 1 CDR3 Curated"],
                "junction_aa": row["Chain 1 junction_aa"],
                "junction": None,
                "consensus_count": None,
                "v_call": row["Chain 1 Curated V Gene"],
                "d_call": None,
                "j_call": row["Chain 1 Curated J Gene"],
                "productive": True,
            }
        )
        chain2.update(
            {
                "locus": receptor_dict[row["Chain 2 Type"]],
                "cdr3_aa": row["Chain 2 CDR3 Curated"],
                "junction_aa": row["Chain 2 junction_aa"],
                "junction": None,
                "consensus_count": None,
                "v_call": row["Chain 2 Curated V Gene"],
                "d_call": row["Chain 2 Curated D Gene"],
                "j_call": row["Chain 2 Curated J Gene"],
                "productive": True,
            }
        )
        for chain_dict in [chain1, chain2]:
            # Since IEDB does not distinguish between lambda and kappa light chains, we need
            # to call them from the gene names
            if chain_dict["locus"] is None:
                chain_dict["locus"] = _infer_locus_from_gene_names(chain_dict, keys=("v_call", "d_call", "j_call"))
            cell.add_chain(chain_dict)

        tcr_cells.append(cell)

    logging.info("Converting to AnnData object")
    iedb = from_airr_cells(tcr_cells)
    iedb_df = iedb_df.set_index(iedb.obs.index)

    iedb.uns["DB"] = {"name": "IEDB", "date_downloaded": datetime.now().isoformat()}
    index_chains(iedb)

    # store cache
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    iedb.write_h5ad(cast(os.PathLike, cache_path))

    return iedb
