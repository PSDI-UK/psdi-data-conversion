import igraph as ig
import matplotlib.pyplot as plt

CONVERTER_OB = "Open Babel"
CONVERTER_ATO = "Atomsk"

# Construct a graph with 4 vertices
l_names = ["MOLDY", "CIF", "PDB", "InChI"]
edges = [(0, 1), (0, 2), (1, 2), (1, 2), (1, 3), (2, 3)]
l_converters = [CONVERTER_ATO, CONVERTER_ATO, CONVERTER_ATO, CONVERTER_OB, CONVERTER_OB, CONVERTER_OB]
g = ig.Graph(len(l_names), edges, vertex_attrs={"name": l_names}, edge_attrs={"converter": l_converters})

# Set title for the graph
g["title"] = "Example conversions"

# Plot in matplotlib
# Note that attributes can be set globally (e.g. vertex_size), or set individually using arrays (e.g. vertex_color)
fig, ax = plt.subplots(figsize=(5, 5))
ig.plot(
    g,
    target=ax,
    # layout="circle",  # print nodes in a circular layout
    vertex_size=30,
    vertex_color="steelblue",
    vertex_frame_width=4.0,
    vertex_frame_color="white",
    vertex_label=g.vs["name"],
    vertex_label_size=7.0,
    # edge_width=[2 if married else 1 for married in g.es["married"]],
    edge_color=["#7142cf" if converter == CONVERTER_ATO else "#AAA" for converter in g.es["converter"]]
)

plt.show()

exit()

# Save the graph as an image file
fig.savefig('social_network.png')
fig.savefig('social_network.jpg')
fig.savefig('social_network.pdf')

# Export and import a graph as a GML file.
g.save("social_network.gml")
g = ig.load("social_network.gml")
