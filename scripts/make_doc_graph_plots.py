import os

import igraph as ig
import matplotlib.pyplot as plt

from psdi_data_conversion.database import get_database

CONVERTER_OB = "Open Babel"
CONVERTER_ATO = "Atomsk"

doc_dir = "./doc"
img_dir = os.path.join(doc_dir, "img")
os.makedirs(img_dir, exist_ok=True)

# Construct a graph with 4 vertices
l_names = ["MOLDY", "CIF", "PDB", "InChI"]
edges = [(0, 1), (0, 2), (1, 2), (1, 2), (1, 3), (2, 3)]
l_converters = [CONVERTER_ATO, CONVERTER_ATO, CONVERTER_ATO, CONVERTER_OB, CONVERTER_OB, CONVERTER_OB]
g = ig.Graph(len(l_names), edges, vertex_attrs={"label": l_names}, edge_attrs={"label": l_converters})

# Set title for the graph
g["title"] = "Example conversions"

# Set up the desired layout of the graph
layout = g.layout(layout="grid")
layout.rotate(-45)

# Plot in matplotlib
fig, ax = plt.subplots(figsize=(5, 5))
ig.plot(
    g,
    target=ax,
    layout=layout,
    vertex_size=30,
    vertex_color="steelblue",
    vertex_frame_width=4.0,
    vertex_frame_color="white",
    vertex_label_size=18,
    vertex_label_dist=1.5,
    edge_width=2,
    edge_color=["red" if converter == CONVERTER_ATO else "blue" if converter == CONVERTER_OB
                else "green" for converter in g.es["label"]]
)

# plt.show()

# Save the graph as an image file
fig.savefig(os.path.join(img_dir, "simple_graph.svg"))

# Add a node and edge
g.add_vertex(label="Custom")
g.add_edge(3, 4, label="Custom")
layout.append((2, 0))

# Plot the new graph
fig, ax = plt.subplots(figsize=(7, 5))
ig.plot(
    g,
    target=ax,
    layout=layout,
    vertex_size=30,
    vertex_color="steelblue",
    vertex_frame_width=4.0,
    vertex_frame_color="white",
    vertex_label_size=18,
    vertex_label_dist=1.5,
    edge_width=2,
    edge_color=["red" if converter == CONVERTER_ATO else "blue" if converter == CONVERTER_OB
                else "green" for converter in g.es["label"]]
)

# Save the new graph
fig.savefig(os.path.join(img_dir, "simple_graph_with_custom.svg"))

# Now let's make a graph of the actual formats and conversions
full_graph = get_database().conversions_table.graph
fig, ax = plt.subplots(figsize=(8, 8))

ig.plot(
    full_graph,
    target=ax,
    layout="auto",
    vertex_size=2,
    vertex_color="black",
    vertex_label_size=4,
    vertex_label_dist=0,
    edge_width=0.1,
    edge_color=["red" if converter == "atomsk" else "blue" if converter == "openbabel"
                else "green" if converter == "c2x" else "grey" for converter in full_graph.es["name"]],
    edge_arrow_size=0,
    edge_arrow_width=0,
)

plt.show()
