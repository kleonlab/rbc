
from pathlib import Path

import anndata as ad
import pandas as pd

DATA_DIR = Path("data/competition_support_set")
GENE_NAMES = DATA_DIR / "gene_names.csv"
TRAINING_H5 = DATA_DIR / "competition_train.h5"
OUTPUT_SAMPLE_H5AD = DATA_DIR / "control_cells_sample.h5ad"
CONTROL_LABEL = "non-targeting"
SAMPLE_SIZE = 100
RANDOM_SEED = 42


def load_gene_names():
    genes = pd.read_csv(
        GENE_NAMES,
        header=None,
        names=["gene_name"],
        dtype=str,
    )
    return genes["gene_name"].tolist()


def sample_control_cells():
    adata = ad.read_h5ad(TRAINING_H5)
    obs = adata.obs

    if "target_gene" not in obs.columns:
        raise KeyError(
            "Column 'target_gene' not found in the AnnData obs table. "
            f"Available columns: {list(obs.columns)}"
        )

    target_series = obs["target_gene"].astype(str).str.lower()
    control_mask = target_series == CONTROL_LABEL
    control_obs = obs.loc[control_mask].copy()

    if control_obs.empty:
        raise ValueError(
            f"No control cells found where target_gene == '{CONTROL_LABEL}'."
        )

    sample_n = min(SAMPLE_SIZE, len(control_obs))
    sampled_obs = control_obs.sample(sample_n, random_state=RANDOM_SEED)
    sampled_adata = adata[sampled_obs.index].copy()
    sampled_adata.write_h5ad(OUTPUT_SAMPLE_H5AD, compression="gzip")

    print(
        f"Wrote {sample_n} control cells (target_gene='{CONTROL_LABEL}') "
        f"to {OUTPUT_SAMPLE_H5AD}"
    )


if __name__ == "__main__":
    gene_list = load_gene_names()
    print(f"Loaded {len(gene_list)} gene names.")
    sample_control_cells()