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

# layout
left, top = 1.0, 3.0     # where the 2x2 grid starts
cell_size = 1.5

# colored cells
for r, c, label, sub, color in cells:
    x = left + c * cell_size
    y = top - (r + 1) * cell_size
    ax.add_patch(Rectangle((x, y), cell_size, cell_size,
                            facecolor=color, edgecolor="black", linewidth=1.5))
    ax.text(x + cell_size / 2, y + cell_size / 2 + 0.15, label,
            ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text(x + cell_size / 2, y + cell_size / 2 - 0.35, sub,
            ha="center", va="center", fontsize=12)


# column labels
ax.text(left + cell_size, top + 0.55, "Actual Label",
        ha="center", va="center", fontsize=14, fontweight="bold")

ax.text(left + cell_size / 2, top + 0.2, "Positive", ha="center", va="center", fontsize=12)
ax.text(left + 1.5 * cell_size, top + 0.2, "Negative", ha="center", va="center", fontsize=12)

# row labels
ax.text(left - 0.85, top - cell_size, "Predicted Label",
        ha="center", va="center", fontsize=14, fontweight="bold", rotation=90)

ax.text(left - 0.2, top - cell_size / 2, "Positive", ha="center", va="center", fontsize=12, rotation=90)
ax.text(left - 0.2, top - 1.5 * cell_size, "Negative", ha="center", va="center", fontsize=12, rotation=90)

# any and all black lines

# outline
ax.add_patch(Rectangle((left, top - 2 * cell_size), 2 * cell_size, 2 * cell_size,
                        fill=False, edgecolor="black", linewidth=2))

# colunm separator
ax.plot([left + cell_size, left + cell_size], [top, top + 0.4], color="black", linewidth=1)
ax.plot([left, left], [top, top + 0.4], color="black", linewidth=1)
ax.plot([left + 2*cell_size -0.01, left + 2*cell_size -0.01], [top, top + 0.4], color="black", linewidth=1)


# row separator
ax.plot([left - 0.4, left], [top - cell_size , top- cell_size], color="black", linewidth=1)
ax.plot([left - 0.4, left], [top , top], color="black", linewidth=1)
ax.plot([left - 0.4, left], [top - 2*cell_size +0.01 , top - 2*cell_size +0.01], color="black", linewidth=1)

# cap row and column
ax.plot([left, left+2*cell_size], [top + 0.4, top + 0.4], color="black", linewidth=1)


# ax.plot([left - 0.55, left - 0.55], [top - 2 * cell_size, top], color="black", linewidth=1.5)

# ax.plot([left, left + 2 * cell_size], [top + 0.3, top + 0.3], color="black", linewidth=1.5)

# ax.plot([left, left + 2 * cell_size], [top - 0.05, top - 0.05], color="black", linewidth=1)

# ax.plot([left - 0.55, left], [top - 0.05, top - 0.05], color="black", linewidth=1)

# ax.plot([left - 0.55, left], [top - cell_size, top - cell_size], color="black", linewidth=1)



plt.tight_layout()
plt.savefig("uhm.png")
# plt.show()
