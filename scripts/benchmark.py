#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser("Summarize the (approximate!) runtime and memory usage of the pipeline from the benchmark files")
parser.add_argument(
    "-o", "--out", default=sys.stdout,
    help="the filename to save the statistics to; defaults to stdout"
)
parser.add_argument(
    "benchmark_directory", nargs='?', type=Path, default="out/benchmark",
    help="the path to the benchmark directory for this pipeline; defaults to 'out/benchmark'"
)
parser.add_argument(
    "workflows", nargs='*',
    help="which workflows should we use?"
)
args = parser.parse_args()

pipeline = Path(args.benchmark_directory)
if not pipeline.is_dir():
    raise ValueError(str(pipeline)+" is not a valid directory")
pipeline = [
    w for w in pipeline.iterdir()
    if not args.workflows or w.stem in args.workflows
]

def bench_rule(rule):
    return pd.DataFrame({
        'max': pd.concat((
            pd.read_csv(f, sep="\t", header=0)
            for f in rule.iterdir()
        )).max()
    }).T

def bench_workflow(workflow):
    return pd.concat({
        rule.stem: bench_rule(rule)
        for rule in workflow.iterdir()
    }, ignore_index=False)

def bench_pipeline(workflows):
    return pd.concat({
        workflow.stem: bench_workflow(workflow)
        for workflow in workflows
    }, ignore_index=False)

pipeline = bench_pipeline(pipeline)[['s', 'max_rss', 'io_out']]
print(pipeline)

print('\ntotals:')
print(pd.Series({'time':pipeline['s'].sum(), 'memory':pipeline['max_rss'].max(), 'disk':pipeline['io_out'].max()}).T)
