import scirpy as ir
import muon as mu
import scanpy as sc
import pandas as pd

def gupta2025():
    adatas_gex = []
    adatas_airr = []
    
    for sample in samples:
        # Load gene expression data
        gex_file = f"data/{sample}/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5"
        adata_gex = sc.read_10x_h5(gex_file)
        adatas_gex.append(adata_gex)
        
        # Load TCR data
        tcr_file = f"data/{sample}/outs/filtered_feature_bc_matrix/filtered_contig_annotations_tcr.csv"
        adata_tcr = ir.io.read_10x_vdj(tcr_file)
        
        # Load BCR data
        bcr_file = f"data/{sample}/outs/filtered_feature_bc_matrix/filtered_contig_annotations_bcr.csv"
        adata_bcr = ir.io.read_10x_vdj(bcr_file)
        
        # Concatenate TCR and BCR data
        adata_airr = adata_tcr.concatenate(adata_bcr, batch_key="receptor_type", batch_categories=["TCR", "BCR"])
        adatas_airr.append(adata_airr)

        # Concatenate all gene expression data
    combined_gex = adatas_gex[0].concatenate(*adatas_gex[1:], batch_key="sample", batch_categories=samples)
    
    # Concatenate all immune-receptor data
    combined_airr = adatas_airr[0].concatenate(*adatas_airr[1:], batch_key="sample", batch_categories=samples)
    
    # Create MuData object
    mdata = mu.MuData({"gex": combined_gex, "airr": combined_airr})
    
    # Use this list of 3k barcodes for consistency with previous versions
    barcodes = pd.read_csv("./barcodes.tsv.csv", header=None)[0].values
    barcodes = pd.Series(barcodes).str.replace("-\\d+$", "", regex=True).values
    
    # Ensure unique observation names
    assert mdata.obs_names.is_unique
    
    # Subset the data to the specified barcodes
    mdata = mdata[barcodes, :].copy()
    
    # Save the subsetted data
    mdata.write_h5mu("snyder2025.h5mu", compression="lzf")
    
    return mdata

# Example usage with multiple samples
samples = ["SH005", "SH006", "SH007"]
mdata = gupta2025(samples)
