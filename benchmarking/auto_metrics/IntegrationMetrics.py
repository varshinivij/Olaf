# --- New metric class using scib-tools ------------------------------------
from scib.metrics import silhouette, silhouette_batch
from scib_metrics.nearest_neighbor import isolated_label_f1 # Corrected import
from typing import Dict

EMBED = "X_scVI"        # change if your embedding key is different
BATCH_KEY = "batch"
LABEL_KEY = "cell_type"

class IntegrationMetric(AutoMetric):
    """
    Compute SCIB integration quality metrics on an AnnData object.
    Returns:
        • batch_silhouette: How well batches mix (lower ≈ better)
        • celltype_silhouette: How well cell types separate (higher ≈ better)
        • isolated_label_f1: Label preservation in isolated clusters (higher ≈ better)
    """

    def metric(self, adata) -> Dict[str, float]:
        # Batch-mixing silhouette (lower is better)
        batch_sil = silhouette_batch(
            adata,
            batch_key=BATCH_KEY,
            embed=EMBED,
            verbose=False
        )

        # Cell-type silhouette (higher is better)
        celltype_sil = silhouette(
            adata,
            label_key=LABEL_KEY,
            embed=EMBED,
            verbose=False
        )

        # Isolated-label F1 (higher is better)
        iso_f1 = isolated_label_f1(
            adata,
            label_key=LABEL_KEY,
            batch_key=BATCH_KEY,
            embed=EMBED,
            verbose=False
        )

        return {
            "Batch Silhouette (↓)": batch_sil,
            "Cell-type Silhouette (↑)": celltype_sil,
            "Isolated-Label F1 (↑)": iso_f1
        }

# --- Example execution -----------------------------------------------------
IntegrationMetric().run(adata)