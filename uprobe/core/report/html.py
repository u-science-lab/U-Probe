"""Public entry point for the universal English HTML report."""
from pathlib import Path
from typing import Optional, Dict, Any
import pandas as pd
from uprobe.core.utils import get_logger
from .compact import save_compact_report

logger = get_logger(__name__)


def save_html_report(df: pd.DataFrame, protocol: Dict[str, Any], output_path: Path,
                     template_type: str = "scientific_report", plot_data=None,
                     csv_filename: Optional[str] = None, raw_df=None) -> Optional[Path]:
    """Legacy template names share one standalone, protocol-driven report.

    summary.report may specify key_metrics and table_columns. Plot data is
    accepted for backward compatibility; distributions use actual final rows.
    """
    try:
        return save_compact_report(df, protocol, output_path, raw_df, csv_filename)
    except Exception:
        logger.exception("Failed to generate HTML report")
        return None
