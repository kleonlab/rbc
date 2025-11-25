from pathlib import Path

import anndata as ad

DATA_DIR = Path("data")
CONTROL_H5AD = DATA_DIR / "competition_support_set/control_cells_sample.h5ad"


from sample_control import load_gene_names


def update_control_targets(goi: str) -> None:
    """Load control subset AnnData and replace target_gene labels."""
    if not CONTROL_H5AD.exists():
        raise FileNotFoundError(
            f"Control AnnData file not found at {CONTROL_H5AD}. "
            "Run scripts/sample_control.py first to generate it."
        )

    adata = ad.read_h5ad(CONTROL_H5AD)

    if "target_gene" not in adata.obs.columns:
        raise KeyError(
            "Column 'target_gene' not found in control AnnData obs table. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    adata.obs["target_gene"] = goi
    adata.write_h5ad(CONTROL_H5AD, compression="gzip")

    print(
        f"Updated {adata.n_obs} control cells in {CONTROL_H5AD} "
        f"with target_gene='{goi}'."
    )


if __name__ == "__main__":
    gene_list = load_gene_names()
    goi = gene_list[1]
    print(gene_list[0:4])
    print(f"Using gene of interest: {goi}")
    update_control_targets(goi)