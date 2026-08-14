import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

metodos = [
    "write_file",
    "_create_alignment"
]

complexidade = [
    50,
    31
]

fig, ax = plt.subplots(figsize=(12, 5))

# Criar barras arredondadas
for i, valor in enumerate(complexidade):
    barra = FancyBboxPatch(
        (i - 0.2, 0),
        0.4,
        valor,
        boxstyle="round,pad=0.02,rounding_size=0.15",
        linewidth=0,
        facecolor="#A98AE8"
    )

    ax.add_patch(barra)

# Limites do gráfico
ax.set_xlim(-0.5, len(metodos) - 0.5)
ax.set_ylim(0, 60)

# Eixo X
ax.set_xticks(range(len(metodos)))
ax.set_xticklabels(metodos)

# Eixo Y
ax.set_yticks([0, 15, 30, 45, 60])

# Linhas horizontais pontilhadas
ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.8,
    alpha=0.3
)

# Grid atrás das barras
ax.set_axisbelow(True)

# Remover bordas
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)

# Remover marcas dos eixos
ax.tick_params(
    axis="both",
    length=0,
    colors="#666666"
)

plt.tight_layout()
plt.show()