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

# O seu código é dinâmico! Basta adicionar o arquivo aqui 
# e ele buscará automaticamente os JSONs correspondentes.
ARQUIVOS = {
    "bigpsl.py": "Bio/Align/bigpsl.py",
    "bigmaf.py": "Bio/Align/bigmaf.py",
    "bed.py": "Bio/Align/bed.py", # <-- ADICIONADO AQUI
}


# ============================================================
# FUNÇÕES
# ============================================================

def carregar_json(nome):
    """Carrega um arquivo JSON da pasta resultados.

    Tenta UTF-8 e UTF-16 porque arquivos gerados pelo
    redirecionamento do PowerShell podem estar em UTF-16.
    """
    caminho = RESULTADOS_DIR / nome

    if not caminho.exists():
        return None

    for encoding in ("utf-8-sig", "utf-16", "utf-8"):
        try:
            with open(caminho, "r", encoding=encoding) as arquivo:
                return json.load(arquivo)
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue

    st.error(f"Não foi possível ler o arquivo JSON: {nome}")
    return None


def classe_radon(complexidade):
    """Converte complexidade para classificação Radon."""
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
    """Executa uma ferramenta externa."""
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
def obter_coesao(arquivo):
    """Executa o Cohesion para o arquivo selecionado."""
    comando = [
        "cohesion",
        "--files",
        arquivo,
        "--verbose",
    ]

    _, saida = executar_comando(comando)

    classes = {}
    classe_atual = None

    for linha in saida.splitlines():
        linha = linha.strip()

        match_classe = re.match(
            r"Class:\s+(.+?)\s+\(\d+:\d+\)",
            linha,
        )

        if match_classe:
            classe_atual = match_classe.group(1)
            continue

        if linha.startswith("Total:") and classe_atual:
            match_total = re.search(
                r"Total:\s+([\d.]+)%",
                linha,
            )

            if match_total:
                classes[classe_atual] = float(
                    match_total.group(1)
                )

                classe_atual = None

    return classes, saida


@st.cache_data
def obter_import_linter():
    """Executa o Import Linter."""
    comando = ["lint-imports"]

    _, saida = executar_comando(comando)

    status = "N/A"
    arquivos = 0
    dependencias = 0

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

    contrato = "BigPsl independence"

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
    """Cria gráfico de barras."""
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


def carregar_resultados_arquivo(nome_arquivo):
    """
    Carrega os resultados de Radon, MI e Pylint
    referentes ao arquivo selecionado.
    """
    # A MÁGICA ESTÁ AQUI: Extrai 'bed' de 'bed.py' automaticamente
    nome = Path(nome_arquivo).stem 

    return {
        "radon_antes": carregar_json(
            f"radon_{nome}_antes.json"
        ),
        "radon_depois": carregar_json(
            f"radon_{nome}_depois.json"
        ),
        "mi_antes": carregar_json(
            f"mi_{nome}_antes.json"
        ),
        "mi_depois": carregar_json(
            f"mi_{nome}_depois.json"
        ),
        "pylint_antes": carregar_json(
            f"pylint_{nome}_antes.json"
        ),
        "pylint_depois": carregar_json(
            f"pylint_{nome}_depois.json"
        ),
    }


def extrair_metodos(radon):
    """Extrai métodos do resultado JSON do Radon (Versão Otimizada)."""
    registros = []

    if not radon:
        return registros

    for blocos in radon.values():
        if not isinstance(blocos, list):
            continue

        for bloco in blocos:
            if bloco.get("type") in ("method", "function"):
                registros.append(
                    {
                        "Método": bloco["name"],
                        "Complexidade": bloco["complexity"],
                        "Classificação": classe_radon(
                            bloco["complexity"]
                        ),
                    }
                )

    return registros


def calcular_media_complexidade(radon):
    """Calcula a complexidade média."""
    metodos = extrair_metodos(radon)

    if not metodos:
        return None

    return sum(
        metodo["Complexidade"]
        for metodo in metodos
    ) / len(metodos)


def obter_mi_valor(dados):
    """Extrai o Maintainability Index do JSON do Radon."""
    if not dados:
        return None

    if isinstance(dados, dict):
        valores = list(dados.values())

        if valores:
            valor = valores[0]

            if isinstance(valor, (int, float)):
                return float(valor)

    return None


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("⚙️ Configurações")

arquivo_nome = st.sidebar.selectbox(
    "Arquivo analisado",
    list(ARQUIVOS.keys()),
)

arquivo = ARQUIVOS[arquivo_nome]

st.sidebar.markdown("---")

st.sidebar.write("**Ferramentas utilizadas:**")
st.sidebar.write("🔵 Radon")
st.sidebar.write("🟢 Pylint")
st.sidebar.write("🟣 Cohesion")
st.sidebar.write("🟠 Import Linter")
st.sidebar.write("📈 Plotly")
st.sidebar.write("🔴 Streamlit")


# ============================================================
# CARREGAMENTO
# ============================================================

resultados = carregar_resultados_arquivo(
    arquivo_nome
)

radon_antes = resultados["radon_antes"]
radon_depois = resultados["radon_depois"]

mi_antes = resultados["mi_antes"]
mi_depois = resultados["mi_depois"]

pylint_antes = resultados["pylint_antes"]
pylint_depois = resultados["pylint_depois"]


# ============================================================
# TÍTULO
# ============================================================

st.title("📊 Dashboard de Qualidade de Software")

st.markdown(
    f"""
    **Projeto:** Análise e refatoração do Biopython

    **Arquivo:** `{arquivo}`
    """
)


# ============================================================
# RADON
# ============================================================

st.subheader("1. Complexidade ciclomática")

media_antes = calcular_media_complexidade(
    radon_antes
)

media_depois = calcular_media_complexidade(
    radon_depois
)

metodos_antes = extrair_metodos(
    radon_antes
)

metodos_depois = extrair_metodos(
    radon_depois
)


# ============================================================
# CARDS
# ============================================================

col1, col2, col3 = st.columns(3)

if metodos_antes:
    maior_antes = max(
        metodos_antes,
        key=lambda x: x["Complexidade"]
    )

    col1.metric(
        "Maior complexidade - Antes",
        maior_antes["Complexidade"],
        maior_antes["Método"],
    )
else:
    col1.metric(
        "Maior complexidade - Antes",
        "N/A",
    )


if metodos_depois:
    maior_depois = max(
        metodos_depois,
        key=lambda x: x["Complexidade"]
    )

    col2.metric(
        "Maior complexidade - Depois",
        maior_depois["Complexidade"],
        maior_depois["Método"],
    )
else:
    col2.metric(
        "Maior complexidade - Depois",
        "N/A",
    )


if media_antes is not None and media_depois is not None:

    reducao = (
        (media_antes - media_depois)
        / media_antes
    ) * 100

    col3.metric(
        "Redução da complexidade média",
        f"{reducao:.2f}%",
    )
else:
    col3.metric(
        "Redução da complexidade média",
        "N/A",
    )


# ============================================================
# COMPARAÇÃO DA COMPLEXIDADE MÉDIA
# ============================================================

if media_antes is not None and media_depois is not None:

    df_media = pd.DataFrame(
        {
            "Estado": [
                "Antes",
                "Depois",
            ],
            "Complexidade": [
                media_antes,
                media_depois,
            ],
        }
    )

    fig_media = criar_grafico_barra(
        x=df_media["Estado"],
        y=df_media["Complexidade"],
        titulo=(
            f"Complexidade média - {arquivo_nome}"
        ),
        nome_eixo_y="Complexidade média",
        formato_texto=".2f",
    )

    st.plotly_chart(
        fig_media,
        use_container_width=True,
    )

else:

    st.warning(
        "Não foram encontrados resultados do Radon "
        "para comparação."
    )


# ============================================================
# MÉTODOS DEPOIS
# ============================================================

st.subheader(
    "2. Complexidade dos métodos após a refatoração"
)

if metodos_depois:

    df_metodos = pd.DataFrame(
        metodos_depois
    )

    df_metodos = df_metodos.sort_values(
        "Complexidade",
        ascending=False,
    )

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

    maior = max(
        df_metodos["Complexidade"]
    )

    fig_metodos.update_layout(
        title=f"Métodos - {arquivo_nome}",
        xaxis=dict(
            title="Complexidade",
            showgrid=True,
            gridcolor="#E5E5E5",
            griddash="dot",
            range=[0, maior * 1.2],
        ),
        yaxis=dict(
            title="",
            autorange="reversed",
            showgrid=False,
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=550,
        showlegend=False,
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


# ============================================================
# COESÃO
# ============================================================

st.subheader("3. Coesão das classes")

coesao, saida_cohesion = obter_coesao(
    arquivo
)

if coesao:

    df_cohesion = pd.DataFrame(
        {
            "Classe": list(coesao.keys()),
            "Coesão (%)": list(coesao.values()),
        }
    )

    fig_cohesion = criar_grafico_barra(
        x=df_cohesion["Classe"],
        y=df_cohesion["Coesão (%)"],
        titulo=f"Coesão - {arquivo_nome}",
        nome_eixo_y="Coesão (%)",
        formato_texto=".2f",
    )

    fig_cohesion.update_yaxes(
        range=[0, 100]
    )

    st.plotly_chart(
        fig_cohesion,
        use_container_width=True,
    )

    st.dataframe(
        df_cohesion,
        use_container_width=True,
        hide_index=True,
    )

else:

    st.warning(
        "Não foi possível obter os resultados "
        "do Cohesion."
    )


# ============================================================
# IMPORT LINTER
# ============================================================

st.subheader("4. Acoplamento arquitetural")

import_linter, saida_linter = (
    obter_import_linter()
)

col1, col2, col3 = st.columns(3)

col1.metric(
    "Contrato",
    import_linter["status"],
)

col2.metric(
    "Arquivos analisados",
    import_linter["arquivos"],
)

col3.metric(
    "Dependências",
    import_linter["dependencias"],
)

if import_linter["status"] == "KEPT":

    st.success(
        "O contrato arquitetural foi mantido."
    )

elif import_linter["status"] == "BROKEN":

    st.error(
        "O Import Linter encontrou uma "
        "violação arquitetural."
    )

else:

    st.warning(
        "Não foi possível interpretar "
        "o resultado do Import Linter."
    )


# ============================================================
# PYLINT
# ============================================================

st.subheader("5. Análise do Pylint")

if pylint_depois:

    df_pylint = pd.DataFrame(
        pylint_depois
    )

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
            titulo=f"Pylint - {arquivo_nome}",
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
        f"Arquivo pylint_{Path(arquivo_nome).stem}_depois.json "
        "não encontrado."
    )


# ============================================================
# RESUMO
# ============================================================

st.subheader("📌 Resumo da análise")

if (
    media_antes is not None
    and media_depois is not None
):

    reducao = (
        (media_antes - media_depois)
        / media_antes
    ) * 100

    st.markdown(
        f"""
        ### `{arquivo_nome}`

        - Complexidade média antes:
          **{media_antes:.2f}**

        - Complexidade média depois:
          **{media_depois:.2f}**

        - Redução da complexidade média:
          **{reducao:.2f}%**

        - Coesão:
          **{len(coesao)} classes analisadas**

        - Acoplamento:
          **{import_linter["status"]}**
        """
    )

else:

    st.info(
        "Ainda não existem dados suficientes "
        "para gerar o resumo completo."
    )