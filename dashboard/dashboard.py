import json

import pandas as pd
import plotly.express as px
import streamlit as st


st.set_page_config(
    page_title="Análise de Qualidade",
    layout="wide"
)

st.title("Análise de Qualidade do Código")
st.write("Resultados do Pylint e Radon")


# =========================
# PYLINT
# =========================

with open("pylint.json", encoding="utf-8") as arquivo:
    pylint_data = json.load(arquivo)

pylint_df = pd.DataFrame(pylint_data)


# =========================
# RADON
# =========================

with open("radon_cc.json", encoding="utf-8") as arquivo:
    radon_data = json.load(arquivo)


# Pegar o primeiro arquivo analisado
arquivo = next(iter(radon_data))

metodos = radon_data[arquivo]

radon_df = pd.DataFrame(metodos)


# =========================
# MÉTRICAS
# =========================

col1, col2, col3 = st.columns(3)

col1.metric(
    "Problemas Pylint",
    len(pylint_df)
)

col2.metric(
    "Métodos analisados",
    len(radon_df)
)

col3.metric(
    "Maior complexidade",
    radon_df["complexity"].max()
)


# =========================
# GRÁFICO RADON
# =========================

st.subheader("Complexidade Ciclomática")

fig = px.bar(
    radon_df,
    x="name",
    y="complexity",
    title="Complexidade dos Métodos"
)

fig.update_layout(
    xaxis_title="Método",
    yaxis_title="Complexidade",
)

st.plotly_chart(
    fig,
    use_container_width=True
)


# =========================
# PYLINT
# =========================

st.subheader("Problemas encontrados pelo Pylint")

contagem = (
    pylint_df["type"]
    .value_counts()
    .reset_index()
)

contagem.columns = ["tipo", "quantidade"]

fig_pylint = px.bar(
    contagem,
    x="tipo",
    y="quantidade",
    title="Problemas por categoria"
)

st.plotly_chart(
    fig_pylint,
    use_container_width=True
)