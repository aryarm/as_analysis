import os
import sys
import pandas as pd


def main(final_path, output_path):
    """
        merge the results.csv.gz files (across all samples) into a single csv
        arguments:
            final_path - the path to the "final" directory
            output_path - the filename of the merged results (optionally gz'd)
    """
    # get the paths to the data
    samples = [os.path.split(x[0])[-1] for x in os.walk(final_path)]
    samples = [sample for sample in samples if sample.startswith("GTEX-")]
    final_paths = []
    remove_samples = set()
    for samp in samples:
        result_path = final_path+"/"+samp+"/result.csv.gz"
        if os.path.isfile(result_path):
            final_paths.append(result_path)
        else:
            remove_samples.add(samp)
    samples = [samp for samp in samples if samp not in remove_samples]

    # import all data into a pandas dataframe
    results = pd.concat(
        (pd.read_csv(res, index_col="gene_name") for res in final_paths),
        keys=samples
    )
    results = results.drop(columns=['gene'])
    results.index.names = ['sample', 'gene']

    if output_path.endswith('.gz'):
        results.to_csv(output_path, compression='gzip')
    else:
        results.to_csv(output_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
