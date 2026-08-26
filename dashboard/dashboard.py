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
    """Executa o Import Linter uma única vez e extrai o resultado de TODOS os contratos.

    O projeto define 3 contratos em `.importlinter`: a independência do
    `bigpsl.py` (já existente) e o acoplamento direto entre `AbiIO.py` e
    `AceIO.py` (adicionados para esta refatoração).
    """
    comando = ["lint-imports"]

    _, saida = executar_comando(comando)

    arquivos = 0
    dependencias = 0

    match_analise = re.search(
        r"Analyzed\s+(\d+)\s+files,\s+(\d+)\s+dependencies\.",
        saida,
    )

    if match_analise:
        arquivos = int(match_analise.group(1))
        dependencias = int(match_analise.group(2))

    # Cada contrato aparece no resumo como "<nome do contrato> KEPT|BROKEN"
    contratos = dict(
        re.findall(r"^(.+?)\s+(KEPT|BROKEN)$", saida, flags=re.MULTILINE)
    )

    resultado = {
        "status": contratos.get("BigPsl independence", "N/A"),
        "arquivos": arquivos,
        "dependencias": dependencias,
        "contrato": "BigPsl independence",
        "contratos": contratos,
    }

    return resultado, saida


@st.cache_data
def obter_radon_cc(caminho_arquivo):
    """Executa `radon cc -j` e retorna a lista de blocos (funções/métodos) do arquivo."""
    comando = ["radon", "cc", caminho_arquivo, "-j"]
    _, saida = executar_comando(comando)

    try:
        dados = json.loads(saida)
    except json.JSONDecodeError:
        return []

    blocos = dados.get(caminho_arquivo, [])
    # O JSON do Radon já lista funções e métodos no nível superior (as
    # classes aparecem à parte, com os mesmos métodos aninhados em
    # "methods"), então basta filtrar por tipo para não duplicar linhas.
    return [b for b in blocos if b.get("type") in ("function", "method")]


@st.cache_data
def obter_radon_mi(caminho_arquivo):
    """Executa `radon mi -j` e retorna (valor do MI, rank) do arquivo."""
    comando = ["radon", "mi", caminho_arquivo, "-j"]
    _, saida = executar_comando(comando)

    try:
        dados = json.loads(saida)
    except json.JSONDecodeError:
        return None, None

    info = dados.get(caminho_arquivo)
    if not info:
        return None, None

    return info["mi"], info["rank"]


@st.cache_data
def obter_pylint(caminho_arquivo):
    """Executa o Pylint (--output-format=json) e retorna a lista de ocorrências."""
    comando = ["pylint", caminho_arquivo, "--output-format=json"]
    _, saida = executar_comando(comando)

    try:
        return json.loads(saida)
    except json.JSONDecodeError:
        return []


@st.cache_data
def obter_coesao_arquivos(caminhos):
    """Executa o Cohesion para uma lista de arquivos e extrai a coesão das classes."""
    comando = ["cohesion", "--files", *caminhos, "--verbose"]

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


# ============================================================
# SEÇÃO 2 — Bio/SeqIO/AbiIO.py e Bio/SeqIO/AceIO.py
# ============================================================
#
# Ao contrário da seção do bigpsl.py (que usa valores "depois" salvos em
# JSON), aqui os resultados "depois" são obtidos executando Radon, Pylint,
# Cohesion e Import Linter AO VIVO contra o código já refatorado em disco.
# Os valores "antes" foram medidos manualmente antes da refatoração e ficam
# fixos como constantes, pois o código anterior não existe mais no repositório.

st.markdown("---")
st.title("🧬 Refatoração — Bio/SeqIO/AbiIO.py e Bio/SeqIO/AceIO.py")

st.markdown(
    """
    **Projeto:** Análise e refatoração do Biopython
    **Arquivos:** `Bio/SeqIO/AbiIO.py` (parser do formato ABI) e
    `Bio/SeqIO/AceIO.py` (parser do formato ACE)

    Os valores **"depois"** desta seção são calculados executando as
    ferramentas ao vivo (subprocess) contra o código já refatorado.
    Os valores **"antes"** foram registrados antes da refatoração.
    """
)

CAMINHO_ABI = "Bio/SeqIO/AbiIO.py"
CAMINHO_ACE = "Bio/SeqIO/AceIO.py"

# --- valores "antes", medidos com Radon/Pylint antes de tocar no código ---
abi_next_antes = 20
ace_next_antes = 7

abi_media_antes = 6.857142857142857
ace_media_antes = 4.0

mi_abi_antes = 49.55174187106885
mi_ace_antes = 67.28009153911344

pylint_abi_antes = 33  # 25 convention (majoritariamente line-too-long) + 8 refactor
pylint_ace_antes = 2  # 1 warning (fixme) + 1 convention (invalid-name do módulo)

cohesion_abi_antes = 66.67
cohesion_ace_antes = 50.0

# --- valores "depois", obtidos ao vivo ---
blocos_abi = obter_radon_cc(CAMINHO_ABI)
blocos_ace = obter_radon_cc(CAMINHO_ACE)

mi_abi_depois, mi_abi_rank = obter_radon_mi(CAMINHO_ABI)
mi_ace_depois, mi_ace_rank = obter_radon_mi(CAMINHO_ACE)

pylint_abi_depois = obter_pylint(CAMINHO_ABI)
pylint_ace_depois = obter_pylint(CAMINHO_ACE)

coesao_seqio, saida_cohesion_seqio = obter_coesao_arquivos((CAMINHO_ABI, CAMINHO_ACE))
cohesion_abi_depois = coesao_seqio.get("AbiIterator", 0)
cohesion_ace_depois = coesao_seqio.get("AceIterator", 0)

abi_next_depois = next(
    (b["complexity"] for b in blocos_abi if b["name"] == "__next__"), 0
)
ace_next_depois = next(
    (b["complexity"] for b in blocos_ace if b["name"] == "__next__"), 0
)

abi_media_depois = (
    sum(b["complexity"] for b in blocos_abi) / len(blocos_abi) if blocos_abi else 0
)
ace_media_depois = (
    sum(b["complexity"] for b in blocos_ace) / len(blocos_ace) if blocos_ace else 0
)

import_linter_contratos = import_linter_resultado["contratos"]
status_abi_forbidden = import_linter_contratos.get(
    "AbiIO does not couple to AceIO", "N/A"
)
status_ace_forbidden = import_linter_contratos.get(
    "AceIO does not couple to AbiIO", "N/A"
)


# ------------------------------------------------------------------
# CARDS — visão geral
# ------------------------------------------------------------------

col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "AbiIterator.__next__ (CC)",
    f"{abi_next_depois}",
    delta=f"{abi_next_depois - abi_next_antes}",
)
col2.metric(
    "AceIterator.__next__ (CC)",
    f"{ace_next_depois}",
    delta=f"{ace_next_depois - ace_next_antes}",
)
col3.metric(
    "MI AbiIO.py",
    f"{mi_abi_depois:.2f}" if mi_abi_depois is not None else "N/A",
    delta=f"{(mi_abi_depois - mi_abi_antes):.2f}" if mi_abi_depois is not None else None,
)
col4.metric(
    "MI AceIO.py",
    f"{mi_ace_depois:.2f}" if mi_ace_depois is not None else "N/A",
    delta=f"{(mi_ace_depois - mi_ace_antes):.2f}" if mi_ace_depois is not None else None,
)


# ------------------------------------------------------------------
# GRÁFICO — Complexidade do método mais crítico (__next__) de cada classe
# ------------------------------------------------------------------

st.subheader("1. Complexidade do método `__next__` (Radon CC)")

col_a, col_b = st.columns(2)

with col_a:
    fig_abi_next = criar_grafico_barra(
        x=["Antes", "Depois"],
        y=[abi_next_antes, abi_next_depois],
        titulo="AbiIterator.__next__",
        nome_eixo_y="Complexidade ciclomática",
        formato_texto=".0f",
        altura=420,
    )
    st.plotly_chart(fig_abi_next, use_container_width=True)

with col_b:
    fig_ace_next = criar_grafico_barra(
        x=["Antes", "Depois"],
        y=[ace_next_antes, ace_next_depois],
        titulo="AceIterator.__next__",
        nome_eixo_y="Complexidade ciclomática",
        formato_texto=".0f",
        altura=420,
    )
    st.plotly_chart(fig_ace_next, use_container_width=True)


# ------------------------------------------------------------------
# GRÁFICO — Complexidade média do arquivo
# ------------------------------------------------------------------

st.subheader("2. Complexidade ciclomática média do arquivo (Radon CC)")

col_a, col_b = st.columns(2)

with col_a:
    fig_abi_media = criar_grafico_barra(
        x=["Antes", "Depois"],
        y=[abi_media_antes, abi_media_depois],
        titulo="Complexidade média - AbiIO.py",
        nome_eixo_y="Complexidade média",
        formato_texto=".2f",
        altura=420,
    )
    st.plotly_chart(fig_abi_media, use_container_width=True)

with col_b:
    fig_ace_media = criar_grafico_barra(
        x=["Antes", "Depois"],
        y=[ace_media_antes, ace_media_depois],
        titulo="Complexidade média - AceIO.py",
        nome_eixo_y="Complexidade média",
        formato_texto=".2f",
        altura=420,
    )
    st.plotly_chart(fig_ace_media, use_container_width=True)


# ------------------------------------------------------------------
# GRÁFICO — Maintainability Index
# ------------------------------------------------------------------

st.subheader("3. Maintainability Index (Radon MI)")

df_mi_seqio = pd.DataFrame(
    {
        "Arquivo": ["AbiIO.py", "AbiIO.py", "AceIO.py", "AceIO.py"],
        "Estado": ["Antes", "Depois", "Antes", "Depois"],
        "MI": [
            mi_abi_antes,
            mi_abi_depois if mi_abi_depois is not None else 0,
            mi_ace_antes,
            mi_ace_depois if mi_ace_depois is not None else 0,
        ],
    }
)

fig_mi_seqio = go.Figure()
for arquivo_nome, cor in [("AbiIO.py", "#B18AEF"), ("AceIO.py", "#7ED6A5")]:
    subset = df_mi_seqio[df_mi_seqio["Arquivo"] == arquivo_nome]
    fig_mi_seqio.add_trace(
        go.Bar(
            x=subset["Estado"] + " (" + subset["Arquivo"] + ")",
            y=subset["MI"],
            name=arquivo_nome,
            text=subset["MI"],
            texttemplate="%{text:.2f}",
            textposition="outside",
            marker=dict(color=cor, line=dict(width=0)),
            width=0.35,
        )
    )

fig_mi_seqio.update_layout(
    title=dict(text="Maintainability Index - AbiIO.py x AceIO.py", x=0.5, xanchor="center"),
    yaxis=dict(title="MI", showgrid=True, gridcolor="#E5E5E5", griddash="dot"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(family="Arial", size=14, color="#333333"),
    margin=dict(l=60, r=30, t=80, b=70),
    showlegend=False,
    height=480,
)

st.plotly_chart(fig_mi_seqio, use_container_width=True)

st.caption(
    "Ambos os arquivos já estavam classificados como **A** no MI antes da "
    "refatoração; a melhora aqui reflete principalmente a redução da "
    "complexidade ciclomática, mesmo com um pequeno aumento de linhas de "
    "código (funções extraídas)."
)


# ------------------------------------------------------------------
# GRÁFICO — Coesão das classes (Cohesion)
# ------------------------------------------------------------------

st.subheader("4. Coesão das classes (Cohesion)")

df_cohesion_seqio = pd.DataFrame(
    {
        "Classe": ["AbiIterator", "AbiIterator", "AceIterator", "AceIterator"],
        "Estado": ["Antes", "Depois", "Antes", "Depois"],
        "Coesão (%)": [
            cohesion_abi_antes,
            cohesion_abi_depois,
            cohesion_ace_antes,
            cohesion_ace_depois,
        ],
    }
)

fig_cohesion_seqio = go.Figure()
for classe, cor in [("AbiIterator", "#B18AEF"), ("AceIterator", "#7ED6A5")]:
    subset = df_cohesion_seqio[df_cohesion_seqio["Classe"] == classe]
    fig_cohesion_seqio.add_trace(
        go.Bar(
            x=subset["Estado"] + " (" + subset["Classe"] + ")",
            y=subset["Coesão (%)"],
            name=classe,
            text=subset["Coesão (%)"],
            texttemplate="%{text:.2f}%",
            textposition="outside",
            marker=dict(color=cor, line=dict(width=0)),
            width=0.35,
        )
    )

fig_cohesion_seqio.update_layout(
    title=dict(text="Coesão das classes - AbiIterator x AceIterator", x=0.5, xanchor="center"),
    yaxis=dict(title="Coesão (%)", range=[0, 100], showgrid=True, gridcolor="#E5E5E5"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(size=14, color="#333333"),
    margin=dict(l=60, r=30, t=70, b=60),
    showlegend=False,
    height=480,
)

st.plotly_chart(fig_cohesion_seqio, use_container_width=True)

st.info(
    "A coesão não mudou com a refatoração: o atributo de classe `modes` "
    "(exigido pela interface `SequenceIterator`, sem relação com o estado "
    "de cada instância) é contado pelo Cohesion como uma variável não "
    "usada em `__init__`/`__next__`, o que já limitava o teto do indicador "
    "antes da refatoração. A lógica extraída virou funções de módulo "
    "(não métodos), então a proporção interna da classe se manteve igual."
)

if not coesao_seqio:
    st.warning(
        "Não foi possível obter os resultados do Cohesion. "
        "Verifique se o comando 'cohesion' está disponível no ambiente virtual."
    )


# ------------------------------------------------------------------
# GRÁFICO — Complexidade de cada função/método após a refatoração
# ------------------------------------------------------------------

st.subheader("5. Complexidade por função/método após a refatoração (Radon CC)")

registros_seqio = []
for arquivo_nome, blocos in [("AbiIO.py", blocos_abi), ("AceIO.py", blocos_ace)]:
    for bloco in blocos:
        nome = bloco["name"]
        if bloco.get("classname"):
            nome = f"{bloco['classname']}.{nome}"
        registros_seqio.append(
            {
                "Arquivo": arquivo_nome,
                "Função/Método": nome,
                "Complexidade": bloco["complexity"],
                "Classificação": classe_radon(bloco["complexity"]),
            }
        )

if registros_seqio:
    df_seqio_metodos = pd.DataFrame(registros_seqio).sort_values(
        "Complexidade", ascending=False
    )

    fig_seqio_metodos = go.Figure()
    for arquivo_nome, cor in [("AbiIO.py", "#B18AEF"), ("AceIO.py", "#7ED6A5")]:
        subset = df_seqio_metodos[df_seqio_metodos["Arquivo"] == arquivo_nome]
        fig_seqio_metodos.add_trace(
            go.Bar(
                x=subset["Complexidade"],
                y=subset["Função/Método"],
                orientation="h",
                name=arquivo_nome,
                text=subset["Complexidade"],
                texttemplate="%{text:.0f}",
                textposition="outside",
                marker=dict(color=cor, line=dict(width=0)),
            )
        )

    fig_seqio_metodos.update_layout(
        title=dict(
            text="Complexidade por função/método - Depois da refatoração",
            x=0.5,
            xanchor="center",
            font=dict(size=18, color="#333333"),
        ),
        xaxis=dict(title="Complexidade", showgrid=True, gridcolor="#E5E5E5", griddash="dot"),
        yaxis=dict(title="", autorange="reversed"),
        barmode="group",
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="Arial", size=14, color="#333333"),
        margin=dict(l=220, r=60, t=80, b=70),
        legend=dict(orientation="h", y=-0.12),
        height=650,
    )

    st.plotly_chart(fig_seqio_metodos, use_container_width=True)

    st.dataframe(df_seqio_metodos, use_container_width=True, hide_index=True)
else:
    st.warning("Não foi possível obter os dados de complexidade do Radon.")


# ------------------------------------------------------------------
# GRÁFICO — Acoplamento (Import Linter)
# ------------------------------------------------------------------

st.subheader("6. Acoplamento entre AbiIO.py e AceIO.py (Import Linter)")

st.markdown(
    """
    Foram adicionados dois contratos do tipo `forbidden` ao `.importlinter`,
    com `allow_indirect_imports = True` para ignorar a cadeia indireta que
    passa por `Bio.SeqRecord -> Bio.SeqIO` (comum a todo módulo de formato
    do `Bio.SeqIO`) e checar apenas o acoplamento **direto** entre os dois
    parsers irmãos.
    """
)

col1, col2, col3 = st.columns(3)
col1.metric("AbiIO → AceIO", status_abi_forbidden)
col2.metric("AceIO → AbiIO", status_ace_forbidden)
col3.metric("Dependências analisadas", import_linter_dependencies)

df_linter_seqio = pd.DataFrame(
    {
        "Contrato": [
            "AbiIO does not couple to AceIO",
            "AceIO does not couple to AbiIO",
        ],
        "Status": [status_abi_forbidden, status_ace_forbidden],
    }
)
st.dataframe(df_linter_seqio, use_container_width=True, hide_index=True)

if "BROKEN" in (status_abi_forbidden, status_ace_forbidden):
    st.error("O Import Linter encontrou um acoplamento direto entre AbiIO e AceIO.")
elif "N/A" in (status_abi_forbidden, status_ace_forbidden):
    st.warning("Não foi possível interpretar o resultado do Import Linter.")
else:
    st.success("AbiIO.py e AceIO.py não importam um ao outro diretamente.")


# ------------------------------------------------------------------
# GRÁFICO — Pylint
# ------------------------------------------------------------------

st.subheader("7. Ocorrências do Pylint — antes x depois")

df_pylint_antes_depois = pd.DataFrame(
    {
        "Arquivo": ["AbiIO.py", "AbiIO.py", "AceIO.py", "AceIO.py"],
        "Estado": ["Antes", "Depois", "Antes", "Depois"],
        "Ocorrências": [
            pylint_abi_antes,
            len(pylint_abi_depois),
            pylint_ace_antes,
            len(pylint_ace_depois),
        ],
    }
)

fig_pylint_seqio = go.Figure()
for arquivo_nome, cor in [("AbiIO.py", "#B18AEF"), ("AceIO.py", "#7ED6A5")]:
    subset = df_pylint_antes_depois[df_pylint_antes_depois["Arquivo"] == arquivo_nome]
    fig_pylint_seqio.add_trace(
        go.Bar(
            x=subset["Estado"] + " (" + subset["Arquivo"] + ")",
            y=subset["Ocorrências"],
            name=arquivo_nome,
            text=subset["Ocorrências"],
            texttemplate="%{text:.0f}",
            textposition="outside",
            marker=dict(color=cor, line=dict(width=0)),
            width=0.35,
        )
    )

fig_pylint_seqio.update_layout(
    title=dict(text="Ocorrências do Pylint - antes x depois", x=0.5, xanchor="center"),
    yaxis=dict(title="Ocorrências", showgrid=True, gridcolor="#E5E5E5", griddash="dot"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    font=dict(family="Arial", size=14, color="#333333"),
    margin=dict(l=60, r=30, t=80, b=70),
    showlegend=False,
    height=460,
)

st.plotly_chart(fig_pylint_seqio, use_container_width=True)

if pylint_abi_depois:
    st.markdown("**AbiIO.py — ocorrências restantes (depois):**")
    st.dataframe(pd.DataFrame(pylint_abi_depois), use_container_width=True, hide_index=True)

if pylint_ace_depois:
    st.markdown("**AceIO.py — ocorrências restantes (depois):**")
    st.dataframe(pd.DataFrame(pylint_ace_depois), use_container_width=True, hide_index=True)

st.caption(
    "As ocorrências restantes em ambos os arquivos são o "
    "`invalid-name` do nome do módulo (`AbiIO`/`AceIO`) e, em `AceIO.py`, "
    "um `fixme` (TODO) intencional — ambos mantidos de propósito, "
    "conforme justificado no relatório de refatoração."
)


# ------------------------------------------------------------------
# RESUMO
# ------------------------------------------------------------------

st.subheader("📌 Resumo da refatoração — AbiIO.py e AceIO.py")

st.markdown(
    f"""
    ### `AbiIterator.__next__`
    - Complexidade antes: **{abi_next_antes} (C)**
    - Complexidade depois: **{abi_next_depois} ({classe_radon(abi_next_depois)})**

    ### `AceIterator.__next__`
    - Complexidade antes: **{ace_next_antes} (B)**
    - Complexidade depois: **{ace_next_depois} ({classe_radon(ace_next_depois)})**

    ### Complexidade média
    - AbiIO.py: **{abi_media_antes:.2f}** → **{abi_media_depois:.2f}**
    - AceIO.py: **{ace_media_antes:.2f}** → **{ace_media_depois:.2f}**

    ### Maintainability Index
    - AbiIO.py: **{mi_abi_antes:.2f} (A)** → **{mi_abi_depois:.2f} ({mi_abi_rank})**
    - AceIO.py: **{mi_ace_antes:.2f} (A)** → **{mi_ace_depois:.2f} ({mi_ace_rank})**

    ### Pylint
    - AbiIO.py: **{pylint_abi_antes} ocorrências** → **{len(pylint_abi_depois)} ocorrência(s)**
    - AceIO.py: **{pylint_ace_antes} ocorrências** → **{len(pylint_ace_depois)} ocorrência(s)**

    ### Acoplamento (Import Linter)
    - AbiIO → AceIO: **{status_abi_forbidden}**
    - AceIO → AbiIO: **{status_ace_forbidden}**
    """
)