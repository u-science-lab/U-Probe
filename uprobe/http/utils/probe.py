from pathlib import Path
from .utils import run_cmd
import typing as t

def uprobe(protocol_yaml: t.Union[str, Path, dict], genomes_yaml: t.Union[str, Path, dict], output_csv: str, 
           raw_results_csv: t.Optional[str] = None, workdir: str = "./"):
    """
    Function to execute uprobe commands with the provided parameters.
    If raw_results_csv is not provided, the function will not generate the raw results CSV file.
    Args:
        protocol_yaml: Path to the protocol YAML file
        genomes_yaml: Path to the genomes YAML file
        output_csv: Path for the output CSV file
        raw_results_csv: Path for the raw results CSV file
        workdir: Working directory, defaults to current directory
    """
    cmd = [
        "uprobe",
        f"--protocol_yaml", protocol_yaml,
        f"--genomes_yaml", genomes_yaml,
        f"--output_csv", output_csv,
        f"--raw_results_csv", raw_results_csv,
        f"--workdir", workdir
    ]
    run_cmd(cmd)
    
