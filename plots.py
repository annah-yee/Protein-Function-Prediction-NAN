import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---- Data / colors ----
green = "#c9e6b0"
red = "#f2b6b0"

cells = [
    # (row, col, label, sublabel, color)
    (0, 0, "True\nPositive", "(TP)", green),
    (0, 1, "False\nPositive", "(FP)", red),
    (1, 0, "False\nNegative", "(FN)", red),
    (1, 1, "True\nNegative", "(TN)", green),
]

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(0, 4)
ax.set_ylim(0, 4)
ax.set_aspect("equal")
ax.axis("off")

# Layout: leave a margin on the left/top for the row/col labels
left, top = 1.0, 3.0     # where the 2x2 grid starts
cell_size = 1.5

# Draw the 4 colored cells
for r, c, label, sub, color in cells:
    x = left + c * cell_size
    y = top - (r + 1) * cell_size
    ax.add_patch(Rectangle((x, y), cell_size, cell_size,
                            facecolor=color, edgecolor="black", linewidth=1.5))
    ax.text(x + cell_size / 2, y + cell_size / 2 + 0.15, label,
            ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text(x + cell_size / 2, y + cell_size / 2 - 0.35, sub,
            ha="center", va="center", fontsize=12)

# Outer border around the whole grid
ax.add_patch(Rectangle((left, top - 2 * cell_size), 2 * cell_size, 2 * cell_size,
                        fill=False, edgecolor="black", linewidth=2))

# --- Column headers: "Actual Label" spanning Positive/Negative ---
ax.text(left + cell_size, top + 0.55, "Actual Label",
        ha="center", va="center", fontsize=14, fontweight="bold")
ax.plot([left, left + 2 * cell_size], [top + 0.3, top + 0.3], color="black", linewidth=1.5)

ax.text(left + cell_size / 2, top + 0.1, "Positive", ha="center", va="center", fontsize=12)
ax.text(left + 1.5 * cell_size, top + 0.1, "Negative", ha="center", va="center", fontsize=12)
ax.plot([left, left + 2 * cell_size], [top - 0.05, top - 0.05], color="black", linewidth=1)
# vertical divider between Positive/Negative headers
ax.plot([left + cell_size, left + cell_size], [top - 0.05, top + 0.3], color="black", linewidth=1)

# --- Row headers: "Predicted Label" spanning Positive/Negative ---
ax.text(left - 0.85, top - cell_size, "Predicted Label",
        ha="center", va="center", fontsize=14, fontweight="bold", rotation=90)
ax.plot([left - 0.55, left - 0.55], [top - 2 * cell_size, top], color="black", linewidth=1.5)

ax.text(left - 0.25, top - cell_size / 2, "Positive", ha="center", va="center", fontsize=12, rotation=90)
ax.text(left - 0.25, top - 1.5 * cell_size, "Negative", ha="center", va="center", fontsize=12, rotation=90)
ax.plot([left - 0.55, left], [top - 0.05, top - 0.05], color="black", linewidth=1)
# horizontal divider between Positive/Negative row headers
ax.plot([left - 0.55, left], [top - cell_size, top - cell_size], color="black", linewidth=1)

# --- Title ---
ax.text(0, 3.85, "(a) Binary confusion matrix", ha="left", va="center",
        fontsize=15, fontweight="bold")

plt.tight_layout()
# plt.savefig("/mnt/user-data/outputs/confusion_matrix.png", dpi=200, bbox_inches="tight")
print("saved")