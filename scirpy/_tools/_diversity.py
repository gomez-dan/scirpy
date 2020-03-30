import numpy as np
from .._util import _is_na
from anndata import AnnData
from typing import Dict
import pandas as pd
from typing import Union


def alpha_diversity(
    adata: AnnData,
    groupby: str,
    *,
    target_col: str = "clonotype",
    inplace: bool = True,
    key_added: Union[None, str] = None
) -> Dict:
    """Computes the alpha diversity of clonotypes within a group. 

    Ignores NaN values. 

    Parameters
    ----------
    adata
        Annotated data matrix 
    groupby 
        Column of `obs` by which the grouping will be performed. 
    target_col
        Column on which to compute the alpha diversity
    inplace
        If True, add a column to `obs`. Otherwise return a DataFrame
        with the alpha diversities. 
    key_added
        Key under which the alpha diversity will be stored if inplace is True. 
        Defaults to `alpha_diversity_{target_col}`. 

    Returns
    -------
    Depending on the value of inplace returns a DataFrame with the alpha diversity
    for each group or adds a column to `adata.obs`. 
    """
    # Could rely on skbio.math if more variants are required.
    def _shannon_entropy(freq):
        np.testing.assert_almost_equal(np.sum(freq), 1)
        return -np.sum(freq * np.log2(freq))

    tcr_obs = adata.obs.loc[~_is_na(adata.obs[target_col]), :]
    clono_counts = (
        tcr_obs.groupby([groupby, target_col]).size().reset_index(name="count")
    )

    diversity = dict()
    for k in tcr_obs[groupby].unique():
        tmp_counts = clono_counts.loc[clono_counts[groupby] == k, "count"].values
        tmp_freqs = tmp_counts / np.sum(tmp_counts)
        diversity[k] = _shannon_entropy(tmp_freqs)

    if inplace:
        key_added = "alpha_diversity_" + target_col if key_added is None else key_added
        adata.obs[key_added] = adata.obs[groupby].map(diversity)
    else:
        return pd.DataFrame().from_dict(diversity, orient="index")
