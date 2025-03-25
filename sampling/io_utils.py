import pandas as pd
from pathlib import Path

def save_results_to_csv(results, output_file, mode='a'):
    output_file = Path(output_file)
    file_exists = output_file.exists()
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, mode=mode, header=not file_exists, index=False)

    