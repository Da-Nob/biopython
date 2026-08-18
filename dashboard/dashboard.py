import json
import re
import subprocess
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# ============================================================
# CONFIGURAÇÃO
# ============================================================

st.set_page_config(
    page_title="Dashboard de Qualidade - Biopython",
    page_icon="📊",
    layout="wide",
)

BASE_DIR = Path(__file__).parent
RESULTADOS_DIR = BASE_DIR / "resultados"


# ============================================================
# FUNÇÕES
# ============================================================

def carregar_json(nome):
    """Carrega um arquivo JSON da pasta resultados."""
    caminho = RESULTADOS_DIR / nome

    if not caminho.exists():
        return None

    with open(caminho, "r", encoding="utf-8-sig") as arquivo:
        return json.load(arquivo)


def classe_radon(complexidade):
    """Converte a complexidade numérica para a classificação Radon."""
    if complexidade <= 5:
        return "A"

    if complexidade <= 10:
        return "B"

    if complexidade <= 20:
        return "C"

    if complexidade <= 30:
        return "D"

    if complexidade <= 40:
        return "E"

    return "F"


def executar_comando(comando):
    """Executa uma ferramenta externa e retorna sua saída."""
    try:
        resultado = subprocess.run(
            comando,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            cwd=BASE_DIR.parent,
        )

        saida = resultado.stdout
        if resultado.stderr:
            saida += "\n" + resultado.stderr

        return resultado.returncode, saida

    except FileNotFoundError as erro:
        return 1, f"Ferramenta não encontrada: {erro}"


@st.cache_data
def obter_coesao():
    """Executa o Cohesion e extrai a coesão das classes."""
    comando = [
        "cohesion",
        "--files",
        "Bio/Align/bigpsl.py",
        "--verbose",
    ]

    _, saida = executar_comando(comando)

    classes = {}

    classe_atual = None

    for linha in saida.splitlines():
        linha = linha.strip()

        match_classe = re.match(r"Class:\s+(.+?)\s+\(\d+:\d+\)", linha)
        if match_classe:
            classe_atual = match_classe.group(1)
            continue

        if linha.startswith("Total:") and classe_atual:
            match_total = re.search(r"Total:\s+([\d.]+)%", linha)
            if match_total:
                classes[classe_atual] = float(match_total.group(1))
                classe_atual = None

    return classes, saida


@st.cache_data
def obter_import_linter():
    """Executa o Import Linter e extrai o resultado do contrato."""
    comando = ["lint-imports"]

    _, saida = executar_comando(comando)

    status = "N/A"
    arquivos = 0
    dependencias = 0
    contrato = "BigPsl independence"

    match_analise = re.search(
        r"Analyzed\s+(\d+)\s+files,\s+(\d+)\s+dependencies\.",
        saida,
    )

    if match_analise:
        arquivos = int(match_analise.group(1))
        dependencias = int(match_analise.group(2))

    if "Contracts: 1 kept, 0 broken." in saida:
        status = "KEPT"
    elif "Contracts: 0 kept, 1 broken." in saida:
        status = "BROKEN"

    return {
        "contrato": contrato,
        "status": status,
        "arquivos": arquivos,
        "dependencias": dependencias,
    }, saida


def criar_grafico_barra(
    x,
    y,
    titulo,
    nome_eixo_y,
    formato_texto=".2f",
    altura=500,
):
    """Cria um gráfico de barras com o padrão visual do dashboard."""

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=x,
            y=y,
            text=y,
            texttemplate=f"%{{text:{formato_texto}}}",
            textposition="outside",
            marker=dict(
                color="#B18AEF",
                line=dict(width=0),
            ),
            width=0.35,
            hovertemplate=(
                "<b>%{x}</b><br>"
                f"{nome_eixo_y}: %{{y:{formato_texto}}}"
                "<extra></extra>"
            ),
        )
    )

    maior_valor = max(y) if len(y) > 0 else 1

    fig.update_layout(
        title=dict(
            text=titulo,
            x=0.5,
            xanchor="center",
            font=dict(
                size=18,
                color="#333333",
            ),
        ),
        xaxis=dict(
            title="",
            showgrid=False,
            zeroline=False,
            tickfont=dict(
                size=13,
                color="#555555",
            ),
        ),
        yaxis=dict(
            title=nome_eixo_y,
            showgrid=True,
            gridcolor="#E5E5E5",
            griddash="dot",
            zeroline=False,
            range=[0, maior_valor * 1.2],
            tickfont=dict(
                size=12,
                color="#777777",
            ),
            title_font=dict(
                size=13,
                color="#555555",
            ),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(
            family="Arial",
            size=14,
            color="#333333",
        ),
        margin=dict(
            l=60,
            r=30,
            t=80,
            b=70,
        ),
        showlegend=False,
        height=altura,
    )

    return fig


# ============================================================
# CARREGAMENTO
# ============================================================

radon_depois = carregar_json("radon_depois.json")
mi_depois = carregar_json("mi_depois.json")
pylint_depois = carregar_json("pylint_depois.json")


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("⚙️ Configurações")

arquivo = st.sidebar.selectbox(
    "Arquivo analisado",
    ["bigpsl.py"],
)

st.sidebar.markdown("---")

st.sidebar.write("**Ferramentas utilizadas:**")
st.sidebar.write("🔵 Radon")
st.sidebar.write("🟢 Pylint")
st.sidebar.write("🟠 Plotly")
st.sidebar.write("🔴 Streamlit")


# ============================================================
# TÍTULO
# ============================================================

st.title("📊 Dashboard de Qualidade de Software")

st.markdown(
    """
    **Projeto:** Análise e refatoração do Biopython  
    **Arquivo:** `Bio/Align/bigpsl.py`
    """
)


# ============================================================
# DADOS DA REFATORAÇÃO
# ============================================================

# Antes da refatoração
write_antes = 50
media_antes = 21.666666666666668
mi_antes = 32.08

# Depois da refatoração
write_depois = 4
media_depois = 7.375
mi_depois_valor = 24.18


# ============================================================
# CARDS
# ============================================================

col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "write_file - Antes",
    "50",
)

col2.metric(
    "write_file - Depois",
    "4",
    delta="-46",
)

col3.metric(
    "Complexidade média",
    "7,38",
    delta="-14,29",
)

col4.metric(
    "Maintainability Index",
    "24,18",
    delta="-7,90",
)


# ============================================================
# GRÁFICO 1 - WRITE_FILE
# ============================================================

st.subheader("1. Complexidade do `AlignmentWriter.write_file`")

df_write = pd.DataFrame(
    {
        "Estado": [
            "Antes",
            "Depois",
        ],
        "Complexidade": [
            write_antes,
            write_depois,
        ],
    }
)

fig_write = criar_grafico_barra(
    x=df_write["Estado"],
    y=df_write["Complexidade"],
    titulo="Complexidade ciclomática - write_file",
    nome_eixo_y="Complexidade",
    formato_texto=".0f",
)

st.plotly_chart(
    fig_write,
    use_container_width=True,
)


# ============================================================
# GRÁFICO 2 - COMPLEXIDADE MÉDIA
# ============================================================

st.subheader("2. Complexidade média do arquivo")

df_media = pd.DataFrame(
    {
        "Estado": [
            "Antes",
            "Depois",
        ],
        "Complexidade média": [
            media_antes,
            media_depois,
        ],
    }
)

fig_media = criar_grafico_barra(
    x=df_media["Estado"],
    y=df_media["Complexidade média"],
    titulo="Complexidade média - bigpsl.py",
    nome_eixo_y="Complexidade média",
    formato_texto=".2f",
)

st.plotly_chart(
    fig_media,
    use_container_width=True,
)


# ============================================================
# GRÁFICO 3 - MAINTAINABILITY INDEX
# ============================================================

st.subheader("3. Maintainability Index")

df_mi = pd.DataFrame(
    {
        "Estado": [
            "Antes",
            "Depois",
        ],
        "Maintainability Index": [
            mi_antes,
            mi_depois_valor,
        ],
    }
)

fig_mi = criar_grafico_barra(
    x=df_mi["Estado"],
    y=df_mi["Maintainability Index"],
    titulo="Maintainability Index - bigpsl.py",
    nome_eixo_y="Maintainability Index",
    formato_texto=".2f",
)

st.plotly_chart(
    fig_mi,
    use_container_width=True,
)

# ============================================================
# COESÃO
# ============================================================

coesao_depois, saida_cohesion = obter_coesao()

cohesion_writer = coesao_depois.get("AlignmentWriter", 0)
cohesion_iterator = coesao_depois.get("AlignmentIterator", 0)


# ============================================================
# ACOPLAMENTO - IMPORT LINTER
# ============================================================

import_linter_resultado, saida_import_linter = obter_import_linter()

import_linter_status = import_linter_resultado["status"]
import_linter_dependencies = import_linter_resultado["dependencias"]
import_linter_files = import_linter_resultado["arquivos"]
import_linter_contract = import_linter_resultado["contrato"]

# ============================================================
# GRÁFICO 4 - COESÃO
# ============================================================

st.subheader("4. Coesão das classes")

df_cohesion = pd.DataFrame(
    {
        "Classe": [
            "AlignmentWriter",
            "AlignmentIterator",
        ],
        "Coesão (%)": [
            cohesion_writer,
            cohesion_iterator,
        ],
    }
)

fig_cohesion = go.Figure()

fig_cohesion.add_trace(
    go.Bar(
        x=df_cohesion["Classe"],
        y=df_cohesion["Coesão (%)"],
        text=df_cohesion["Coesão (%)"],
        texttemplate="%{text:.2f}%",
        textposition="outside",
        marker=dict(
            color="#B18AEF",
            line=dict(width=0),
        ),
        hovertemplate=(
            "<b>%{x}</b><br>"
            "Coesão: %{y:.2f}%"
            "<extra></extra>"
        ),
    )
)

fig_cohesion.update_layout(
    title="Coesão das classes - bigpsl.py",
    xaxis=dict(
        title="",
        showgrid=False,
        zeroline=False,
    ),
    yaxis=dict(
        title="Coesão (%)",
        showgrid=True,
        gridcolor="#E5E5E5",
        zeroline=False,
        range=[0, 100],
    ),
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(
        size=14,
        color="#333333",
    ),
    margin=dict(
        l=60,
        r=30,
        t=70,
        b=60,
    ),
    showlegend=False,
)

fig_cohesion.update_traces(
    width=0.35,
)

st.plotly_chart(
    fig_cohesion,
    use_container_width=True,
)

if not coesao_depois:
    st.warning(
        "Não foi possível obter os resultados do Cohesion. "
        "Verifique se o comando 'cohesion' está disponível no ambiente virtual."
    )


# ============================================================
# GRÁFICO 5 - MÉTODOS DEPOIS DA REFATORAÇÃO
# ============================================================

st.subheader(
    "5. Complexidade dos métodos após a refatoração"
)

if radon_depois:

    registros = []

    for arquivo_nome, metodos in radon_depois.items():

        for metodo in metodos:

            if metodo.get("type") == "method":

                registros.append(
                    {
                        "Método": metodo["name"],
                        "Complexidade": metodo["complexity"],
                        "Classificação": classe_radon(
                            metodo["complexity"]
                        ),
                    }
                )

    if registros:

        df_metodos = pd.DataFrame(registros)

        df_metodos = df_metodos.sort_values(
            "Complexidade",
            ascending=False,
        )

        # Gráfico horizontal para facilitar a leitura dos nomes dos métodos
        fig_metodos = go.Figure()

        fig_metodos.add_trace(
            go.Bar(
                x=df_metodos["Complexidade"],
                y=df_metodos["Método"],
                orientation="h",
                text=df_metodos["Complexidade"],
                texttemplate="%{text:.0f}",
                textposition="outside",
                marker=dict(
                    color="#B18AEF",
                    line=dict(width=0),
                ),
                width=0.35,
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "Complexidade: %{x:.0f}"
                    "<extra></extra>"
                ),
            )
        )

        maior_complexidade = (
            max(df_metodos["Complexidade"])
            if not df_metodos.empty
            else 1
        )

        fig_metodos.update_layout(
            title=dict(
                text="Complexidade dos métodos - Depois",
                x=0.5,
                xanchor="center",
                font=dict(size=18, color="#333333"),
            ),
            xaxis=dict(
                title="Complexidade",
                showgrid=True,
                gridcolor="#E5E5E5",
                griddash="dot",
                zeroline=False,
                range=[0, maior_complexidade * 1.2],
                tickfont=dict(size=12, color="#777777"),
                title_font=dict(size=13, color="#555555"),
            ),
            yaxis=dict(
                title="",
                showgrid=False,
                zeroline=False,
                tickfont=dict(size=13, color="#555555"),
                autorange="reversed",
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            font=dict(family="Arial", size=14, color="#333333"),
            margin=dict(l=180, r=60, t=80, b=70),
            showlegend=False,
            height=550,
        )

        st.plotly_chart(
            fig_metodos,
            use_container_width=True,
        )

        st.dataframe(
            df_metodos,
            use_container_width=True,
            hide_index=True,
        )

    else:

        st.warning(
            "Não foram encontrados métodos no JSON do Radon."
        )

else:

    st.warning(
        "O arquivo radon_depois.json não foi encontrado."
    )


# ============================================================
# GRÁFICO 6 - IMPORT LINTER - ACOPLAMENTO
# ============================================================

st.subheader("6. Acoplamento arquitetural")

col1, col2, col3 = st.columns(3)

col1.metric(
    "Status do contrato",
    import_linter_status,
)

col2.metric(
    "Arquivos analisados",
    import_linter_files,
)

col3.metric(
    "Dependências analisadas",
    import_linter_dependencies,
)

st.markdown(
    """
    **Contrato:** `{import_linter_contract}`

    O contrato foi classificado como **KEPT**, indicando que
    a regra arquitetural definida no Import Linter foi respeitada.
    """
)

df_linter = pd.DataFrame(
    {
        "Contrato": ["BigPsl independence"],
        "Status": [import_linter_status],
        "Arquivos analisados": [import_linter_files],
        "Dependências": [import_linter_dependencies],
    }
)

st.dataframe(
    df_linter,
    use_container_width=True,
    hide_index=True,
)

if import_linter_status == "BROKEN":
    st.error(
        "O Import Linter encontrou uma violação no contrato arquitetural."
    )
elif import_linter_status == "N/A":
    st.warning(
        "Não foi possível interpretar o resultado do Import Linter."
    )
else:
    st.success(
        "O contrato arquitetural foi mantido pelo Import Linter."
    )

# ============================================================
# GRÁFICO 7 - PYLINT
# ============================================================

st.subheader("7. Análise do Pylint")

if pylint_depois is not None:

    df_pylint = pd.DataFrame(pylint_depois)

    if not df_pylint.empty:

        contagem = (
            df_pylint["type"]
            .value_counts()
            .reset_index()
        )

        contagem.columns = [
            "Tipo",
            "Quantidade",
        ]

        fig_pylint = criar_grafico_barra(
            x=contagem["Tipo"],
            y=contagem["Quantidade"],
            titulo="Problemas encontrados pelo Pylint",
            nome_eixo_y="Quantidade",
            formato_texto=".0f",
        )

        st.plotly_chart(
            fig_pylint,
            use_container_width=True,
        )

        st.dataframe(
            df_pylint,
            use_container_width=True,
            hide_index=True,
        )

    else:

        st.success(
            "O Pylint não encontrou problemas."
        )

else:

    st.warning(
        "O arquivo pylint_depois.json não foi encontrado."
    )


# ============================================================
# RESUMO
# ============================================================

st.subheader("📌 Resumo da refatoração")

st.markdown(
    f"""
    ### `AlignmentWriter.write_file`

    - Complexidade antes: **50 (F)**
    - Complexidade depois: **4 (A)**
    - Redução da complexidade: **92%**

    ### Complexidade média

    - Antes: **21,67 (D)**
    - Depois: **7,38 (B)**
    - Redução: **65,96%**

    ### Maintainability Index

    - Antes: **32,08 (A)**
    - Depois: **24,18 (A)**

    O índice permaneceu na classificação **A**, embora seu
    valor numérico tenha diminuído.
    """
)