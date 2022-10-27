"""Calculate genealogical nearest neighbours from tree sequence.

Calculate genealogical nearest neighbours (GNN) from tree sequence TS
file. Calculations can be based on either 'population' (default) or
'individual' mode, which defines at what level samples are grouped.

"""
import json

import click
import pandas as pd
import tskit
from tswf.cli import pass_environment


def make_unique(x):
    indices = []
    seen = {}
    for y in x:
        if y not in seen:
            seen[y] = 0
        else:
            seen[y] = seen[y] + 1
        indices.append(seen[y])
    return indices


@click.command(help=__doc__)
@click.argument("ts", type=click.Path(exists=True))
@click.option(
    "--mode", type=click.Choice(["individual", "population"]), default="individual"
)
@click.option("--output-file", type=click.Path())
@pass_environment
def main(env, ts, mode, output_file):
    if output_file is None:
        output_file = f"{ts}.gnn.csv"

    ts = tskit.load(ts)

    if mode == "population":
        samples_listed_by_group = [
            ts.samples(population=pop_id) for pop_id in range(ts.num_populations)
        ]
    elif mode == "individual":
        samples_listed_by_group = [ind.nodes for ind in ts.individuals()]

    gnn = ts.genealogical_nearest_neighbours(ts.samples(), samples_listed_by_group)

    sample_nodes = [ts.node(n) for n in ts.samples()]
    sample_node_ids = [n.id for n in sample_nodes]
    sample_names = [
        json.loads(ts.individual(n.individual).metadata)["SM"] for n in sample_nodes
    ]

    sample_names_unique = list(
        map(lambda x: f"{x[0]}/{x[1]}", zip(sample_names, make_unique(sample_names)))
    )
    sample_node_pop_ids = [ts.population(n.population).id for n in sample_nodes]

    sample_node_pops = [
        json.loads(ts.population(i).metadata.decode())["population"]
        for i in sample_node_pop_ids
    ]

    if mode == "population":
        columns = [json.loads(p.metadata)["population"] for p in ts.populations()]
    elif mode == "individual":
        columns = [json.loads(ind.metadata)["SM"] for ind in ts.individuals()]

    # Could possibly add metadata information in populations here
    gnn_table = pd.DataFrame(
        data=gnn,
        index=[
            pd.Index(sample_node_ids, name="sample_node_id"),
            pd.Index(sample_names, name="sample_name"),
            pd.Index(sample_names_unique, name="sample_name_unique"),
            pd.Index(sample_node_pop_ids, name="sample_node_population_id"),
            pd.Index(sample_node_pops, name="sample_node_population"),
        ],
        columns=columns,
    )

    # Save gnn table
    gnn_table.to_csv(output_file, header=True)  # noqa: F821
