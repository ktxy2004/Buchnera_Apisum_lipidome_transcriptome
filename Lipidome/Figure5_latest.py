import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.patches as mpatches

# --------------------------------------------------
# 1 Define pathway edges
# --------------------------------------------------

edges = [
    ("DHAP", "G3P"),
    ("G3P", "LPA"),
    ("LPA", "PA"),
    ("PA", "CDP-DAG"),
    ("CDP-DAG", "PI"),
    ("CDP-DAG", "PG"),
    ("CDP-DAG", "PS"),
    ("PG", "CL"),
    ("PA", "DG"),
    ("DG", "PC"),
    ("DG", "PE"),
    ("PC", "LPC"),
    ("PE", "LPE")
]

G = nx.DiGraph()
G.add_edges_from(edges)

# --------------------------------------------------
# 2 Node positions
# --------------------------------------------------

pos = {
    "DHAP": (0, 3),
    "G3P": (1.6, 3),
    "LPA": (3.2, 3),
    "PA": (4.8, 3),
    "CDP-DAG": (5.1, 4.2),
    "PI": (7.3, 4.9),
    "PG": (7.3, 4.2),
    "PS": (7.3, 3.6),
    "CL": (9, 4.2),
    "DG": (6.3, 2),
    "PC": (7.6, 2.6),
    "LPC": (9, 2.6),
    "PE": (7.6, 1.5),
    "LPE": (9, 1.5)
}

# --------------------------------------------------
# 3 Node colors
# --------------------------------------------------

host_nodes = ["DHAP", "G3P", "LPA", "PA", "CDP-DAG"]

buchnera_enriched_nodes = [
    "DG", "PI", "PG", "CL", "PC", "PE", "PS", "LPC", "LPE"
]

node_colors = [
    "white" if n in host_nodes else "#FFEFA1"
    for n in G.nodes()
]

# --------------------------------------------------
# 4 Draw figure
# --------------------------------------------------

plt.figure(figsize=(10, 6))

# Light gray background = bacteriocyte context
plt.axvspan(-0.5, 10, color="#f2f2f2", alpha=0.3)

nx.draw_networkx_edges(
    G, pos,
    width=2,
    arrows=True,
    edge_color="#555555",
    arrowsize=14
)

nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    edgecolors="black",
    node_size=2200
)

nx.draw_networkx_labels(
    G, pos,
    font_size=10,
    font_family="DejaVu Sans",
    font_color="black"
)

# --------------------------------------------------
# 5 Enzyme labels
# --------------------------------------------------

enzyme_positions = {
    ("DHAP", "G3P"): ("GPDH ↑", "red", 0.25),
    ("G3P", "LPA"): ("GPAT1 ↑", "red", 0.25),
    ("PA", "CDP-DAG"): ("CDS ↑", "red", 0.10),
    ("PA", "DG"): ("PPAP2A ↑", "red", -0.30),
    ("PG", "CL"): ("cls (Buchnera)", "black", 0.25),
    ("PC", "LPC"): ("GIIIspla2 ↓", "blue", 0.25),
    ("PE", "LPE"): ("GIIIspla2 ↓", "blue", 0.25)
}

for (start, end), (label, color, offset) in enzyme_positions.items():
    x1, y1 = pos[start]
    x2, y2 = pos[end]

    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2

    plt.text(
        mid_x,
        mid_y + offset,
        label,
        fontsize=10,
        color=color,
        ha="center",
        fontweight="bold",
        fontstyle="italic" if "cls" in label else "normal"
    )

# Fine-tune cls placement
for text in plt.gca().texts:
    if "cls" in text.get_text():
        x, y = text.get_position()
        text.set_position((x + 0.4, y + 0.15))

# Fine-tune CDS label
for text in plt.gca().texts:
    if text.get_text().strip() == "CDS ↑":
        x, y = text.get_position()
        text.set_position((x - 0.25, y))

# --------------------------------------------------
# 6 Titles and legend
# --------------------------------------------------

plt.text(
    2.5, 5.2,
    "Host transcriptional regulation of phospholipid metabolism",
    fontsize=12,
    ha="center",
    color="#444444",
    fontweight="bold"
)

# Yellow nodes now explicitly mean Buchnera enrichment
buchnera_patch = mpatches.Patch(
    facecolor="#FFEFA1",
    edgecolor="black",
    label="Buchnera-enriched lipid species"
)

plt.legend(
    handles=[buchnera_patch],
    loc="lower left",
    frameon=False,
    fontsize=9,
    handlelength=1.3,
    labelspacing=0.7,
    borderaxespad=0.6,
    prop={"weight": "bold"}
)

# --------------------------------------------------
# 7 Export
# --------------------------------------------------

plt.axis("off")
plt.tight_layout()

plt.savefig(
    "aphid_buchnera_pathway_final.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "aphid_buchnera_pathway_final_transparent.svg",
    dpi=300,
    transparent=True,
    bbox_inches="tight"
)

plt.show()
