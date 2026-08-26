# Relatório de Refatoração — `Bio/SeqIO/AbiIO.py` e `Bio/SeqIO/AceIO.py`

**Projeto:** Biopython — disciplina de Teste e Qualidade de Software
**Escopo:** refatoração dos parsers `AbiIO.py` (formato ABI) e `AceIO.py` (formato ACE) do `Bio.SeqIO`
**Ferramentas:** Radon (complexidade ciclomática e Maintainability Index), `cognitive-complexity`
(complexidade cognitiva), Pylint (code smells e duplicação), Import Linter (acoplamento), Cohesion (coesão)
**Ambiente:** todas as ferramentas foram executadas via `venv/Scripts/*.exe`, a partir da raiz do
repositório `biopython/`, contra os arquivos reais do projeto (nenhum teste automatizado foi criado —
apenas refatoração, conforme solicitado).

> **Nota sobre "complexidade cognitiva":** o Radon **não calcula** complexidade cognitiva — ele calcula
> apenas complexidade ciclomática (`radon cc`) e Maintainability Index (`radon mi`). Para medir a
> complexidade cognitiva (métrica da SonarSource, que penaliza aninhamento/`elif`/`break` de forma
> diferente da ciclomática) foi instalado o pacote `cognitive-complexity`, que percorre a AST de cada
> função e devolve um score. Isso é explicado com mais detalhe na seção 2.2.

---

## 1. Placar por arquivo — cada métrica, antes e depois

Visão rápida: para cada arquivo, o que cada uma das 7 métricas avaliadas marcava antes da refatoração e o
que passou a marcar depois. Os números "do arquivo" (CC, cognitiva) são a **média entre todas as
funções/métodos do arquivo**; ao lado de cada um vai, entre parênteses, o método que mais pesava nessa
média, para quem quiser o detalhe. A tabela método a método completa está nas seções 3 e 4.

### `Bio/SeqIO/AbiIO.py`

| Métrica | Ferramenta | Antes | Depois | Método que mais pesava |
|---|---|---|---|---|
| Complexidade ciclomática (CC) | Radon (`radon cc`) | **6,86** (média, rank B) | **4,09** (média, rank A) | `AbiIterator.__next__`: 20 (C) → 5 (A) |
| Complexidade cognitiva | `cognitive-complexity` | **7,86** (média) | **3,09** (média) | `AbiIterator.__next__`: 20 → 5 |
| Índice de manutenibilidade (MI) | Radon (`radon mi`) | **49,55** (rank A) | **52,21** (rank A) | — (métrica é por arquivo inteiro) |
| Code smells | Pylint | **33 ocorrências** (nota 7,92/10) | **1 ocorrência** (nota 9,93/10) | `__next__` concentrava 6 dessas 33 (too-many-branches/statements/locals) |
| Duplicação de código | Pylint (`R0801`) + leitura manual | **1 bloco duplicado** (`try/except AttributeError` do nome do arquivo, repetido 2× dentro de `__next__`) — Pylint não acusa por ser curto e intra-arquivo | **0 blocos duplicados** | trecho extraído para a função `_abi_file_stem` |
| Coesão | Cohesion | **66,67%** (classe `AbiIterator`) | **66,67%** (sem mudança) | limitação do atributo `modes`, explicada no §2.5 |
| Acoplamento | Import Linter | **sem contrato definido** | **KEPT** (não importa `AceIO.py` diretamente) | contrato `AbiIO does not couple to AceIO` |

### `Bio/SeqIO/AceIO.py`

| Métrica | Ferramenta | Antes | Depois | Método que mais pesava |
|---|---|---|---|---|
| Complexidade ciclomática (CC) | Radon (`radon cc`) | **4,00** (média, rank A) | **2,50** (média, rank A) | `AceIterator.__next__`: 7 (B) → 2 (A) |
| Complexidade cognitiva | `cognitive-complexity` | **3,50** (média) | **1,50** (média) | `AceIterator.__next__`: 7 → 1 |
| Índice de manutenibilidade (MI) | Radon (`radon mi`) | **67,28** (rank A) | **67,81** (rank A) | — (métrica é por arquivo inteiro) |
| Code smells | Pylint | **2 ocorrências** (nota 9,41/10) | **2 ocorrências** (nota 9,46/10) | nenhuma — as 2 ocorrências são intencionais (ver §2.3) |
| Duplicação de código | Pylint (`R0801`) + leitura manual | **0 blocos duplicados** | **0 blocos duplicados** | arquivo já era pequeno e sem repetição |
| Coesão | Cohesion | **50,00%** (classe `AceIterator`) | **50,00%** (sem mudança) | limitação do atributo `modes`, explicada no §2.5 |
| Acoplamento | Import Linter | **sem contrato definido** | **KEPT** (não importa `AbiIO.py` diretamente) | contrato `AceIO does not couple to AbiIO` |

**Como ler isso:** em `AbiIO.py`, praticamente toda a melhora (CC, cognitiva, code smells) vem de um único
método, `__next__` — era ele que arrastava a média do arquivo para cima. Em `AceIO.py` a história é a
mesma em escala menor: um arquivo já simples, com `__next__` sendo o único ponto realmente complexo, e que
ficou ainda mais simples depois. As duas métricas que **não mudam** (coesão e as 2 ocorrências restantes de
Pylint em `AceIO.py`) não mudam por decisão consciente, não por limitação da refatoração — o motivo de cada
uma está explicado nas seções 2.5 e 2.3, e resumido na tabela acima.

---

## 2. Como cada ferramenta foi usada e o que ela retornou

### 2.1 Radon — Complexidade Ciclomática (CC) e Maintainability Index (MI)

**Como foi usado:**
```bash
radon cc Bio/SeqIO/AbiIO.py Bio/SeqIO/AceIO.py -s      # complexidade ciclomática, por função/classe
radon mi Bio/SeqIO/AbiIO.py Bio/SeqIO/AceIO.py -s      # maintainability index, por arquivo
radon cc Bio/SeqIO/AbiIO.py Bio/SeqIO/AceIO.py -j      # saída em JSON (usada pelo dashboard)
```

**Como a métrica funciona:** o Radon percorre a AST e conta 1 ponto de partida + 1 para cada estrutura
que cria um novo caminho de execução (`if`, `elif`, `for`, `while`, `and`/`or` em condições, `except`,
`with`, comprehensions, etc). O resultado é classificado em letras:

| Rank | Faixa de CC | Interpretação |
|---|---|---|
| A | 1–5  | simples, baixo risco |
| B | 6–10 | pouco complexo |
| C | 11–20 | complexo, moderadamente arriscado |
| D | 21–30 | complexo, arriscado |
| E/F | 31+ | muito complexo, alto risco |

O **Maintainability Index (MI)** combina CC, Halstead Volume e linhas de código numa fórmula única,
normalizada de 0 a 100 (quanto maior, mais fácil de manter). O Radon classifica em A (MI ≥ 20 — fácil
manutenção), B (10–19) e C (< 10 — difícil manutenção).

**Resultado por classe/arquivo:**

| Alvo | CC antes | CC depois | MI antes | MI depois |
|---|---|---|---|---|
| `AbiIterator.__next__` | 20 (C) | 5 (A) | — | — |
| `AbiIterator` (classe, agregado) | 13 (C) | 6 (B) | — | — |
| `_parse_tag_data` (função módulo) | 12 (C) | 5 (A) | — | — |
| `_abi_trim` (função módulo) | 6 (B) | 6 (B) | — | — |
| `Bio/SeqIO/AbiIO.py` (arquivo) | média 6,86 | média 4,09 | 49,55 (A) | 52,21 (A) |
| `AceIterator.__next__` | 7 (B) | 2 (A) | — | — |
| `AceIterator` (classe, agregado) | 5 (A) | 3 (A) | — | — |
| `Bio/SeqIO/AceIO.py` (arquivo) | média 4,00 | média 2,50 | 67,28 (A) | 67,81 (A) |

O MI de `AbiIO.py` subiu de 49,55 para 52,21 principalmente porque a complexidade ciclomática total caiu
bastante; o pequeno aumento de linhas de código (+46, por causa das novas funções extraídas) pesa contra,
mas o efeito da complexidade domina. `AceIO.py` já estava com MI alto (67,28) e subiu pouco (67,81), pois
sua complexidade original já era baixa.

### 2.2 `cognitive-complexity` — Complexidade Cognitiva

**Como foi usado:** script Python que percorre a AST de cada arquivo com `ast.parse` e chama
`get_cognitive_complexity(node)` (da lib `cognitive-complexity`) para cada `FunctionDef`.

**Como a métrica funciona:** diferente da ciclomática (que conta caminhos), a cognitiva **penaliza
aninhamento** — um `if` dentro de outro `if` custa mais do que dois `if` soltos — e ignora estruturas que
não dificultam a leitura (como `and`/`or` puros contam pouco, `switch`/dispatch por dicionário não conta
nada). É pensada para medir o quão difícil é *ler e entender* a função, não só testá-la.

**Resultado por função:**

| Função | Cognitiva antes | Cognitiva depois |
|---|---|---|
| `AbiIterator.__next__` | 20 | 5 |
| `_parse_tag_data` | 16 | 5 |
| `_abi_trim` | 12 | 6 |
| `_abi_parse_header` | 3 | 3 (sem mudança de fluxo, só limpeza de variáveis) |
| `AceIterator.__next__` | 7 | 1 |

A complexidade cognitiva caiu proporcionalmente mais que a ciclomática em quase todos os casos
(`_parse_tag_data`: CC 12→5 mas cognitiva 16→5) porque o dicionário de despacho e a remoção dos `else`
aninhados eliminam justamente o tipo de estrutura que essa métrica penaliza mais (aninhamento e cadeias
`if/elif/else`).

### 2.3 Pylint — Code Smells e Duplicação

**Como foi usado:**
```bash
pylint Bio/SeqIO/AbiIO.py Bio/SeqIO/AceIO.py --output-format=json:resultado.json,text
```
O dashboard executa o mesmo comando ao vivo, por arquivo, com `--output-format=json`.

**Como a métrica funciona:** o Pylint aplica um conjunto de *checkers* estáticos e classifica cada
ocorrência por tipo — `convention` (estilo, ex.: `line-too-long`, `invalid-name`), `refactor` (estrutura,
ex.: `too-many-branches`, `too-many-locals`, `too-many-return-statements`, `no-else-return`) e `warning`
(possíveis bugs, ex.: `fixme`). Ao final soma uma nota de 0 a 10 (10 − penalidades/instrução). Duplicação de
código é o checker `R0801 duplicate-code` (comparação de blocos ≥ `min-similarity-lines` entre arquivos);
não foi disparado entre esses dois arquivos porque eles não compartilham blocos de 4+ linhas idênticas —
mas havia duplicação **dentro** do próprio `AbiIO.py` (o bloco `try/except AttributeError` para extrair o
nome do arquivo, repetido em dois pontos de `__next__`), que o Pylint não sinaliza por padrão (duplicação
intra-arquivo curta), mas que a leitura manual do código identificou e que foi removida na refatoração.

**Resultado por arquivo (contagem de ocorrências por tipo):**

| Arquivo | Antes (total) | Antes por tipo | Depois (total) | Depois por tipo |
|---|---|---|---|---|
| `AbiIO.py` | 33 | `line-too-long`×23, `invalid-name`×2, `too-many-locals`×2, `no-else-return`×2, `no-else-raise`×1, `too-many-branches`×1, `too-many-statements`×1, `too-many-return-statements`×1 | 1 | `invalid-name`×1 (módulo, justificado — ver §3) |
| `AceIO.py` | 2 | `fixme`×1, `invalid-name`×1 | 2 | `fixme`×1, `invalid-name`×1 (ambos mantidos de propósito — ver §4) |

Nota geral do Pylint: `AbiIO.py` foi de **7,92/10** para **9,93/10**; `AceIO.py` já estava em **9,41/10** e
foi para **9,46/10** (não tinha praticamente nada a corrigir).

### 2.4 Import Linter — Acoplamento

**Como foi usado:** o projeto já tinha um contrato de exemplo (`bigpsl-independence`) em `.importlinter`.
Foram adicionados dois contratos novos do tipo `forbidden`, um para cada sentido, com
`allow_indirect_imports = True`:

```ini
[importlinter:contract:abiio-forbidden]
name = AbiIO does not couple to AceIO
type = forbidden
source_modules = Bio.SeqIO.AbiIO
forbidden_modules = Bio.SeqIO.AceIO
allow_indirect_imports = True

[importlinter:contract:aceio-forbidden]
name = AceIO does not couple to AbiIO
type = forbidden
source_modules = Bio.SeqIO.AceIO
forbidden_modules = Bio.SeqIO.AbiIO
allow_indirect_imports = True
```

**Como a métrica funciona:** o Import Linter constrói o grafo de imports de todo o `Bio` e verifica se
existe algum caminho do(s) `source_modules` até os `forbidden_modules`. Com `allow_indirect_imports = True`,
apenas imports **diretos** quebram o contrato — cadeias indiretas (como `AbiIO -> Bio.SeqRecord ->
Bio.SeqIO -> AceIO`, que existe porque `Bio.SeqRecord` importa `Bio.SeqIO` para type hints, e `Bio.SeqIO`
importa todos os parsers de formato) são ignoradas. Sem essa opção, o contrato dava `BROKEN` só por causa
dessa cadeia indireta estrutural do pacote — um falso positivo em relação ao que realmente queríamos medir
(se os dois parsers "irmãos" se importam um ao outro diretamente).

**Resultado:**

```
Analyzed 293 files, 815 dependencies.
BigPsl independence            KEPT
AbiIO does not couple to AceIO KEPT
AceIO does not couple to AbiIO KEPT
Contracts: 3 kept, 0 broken.
```

Os dois parsers não importam um ao outro nem antes nem depois da refatoração — o que era esperado, já que
a refatoração ficou contida dentro de cada arquivo. O contrato serve como guarda continuada: qualquer PR
futuro que crie esse acoplamento vai quebrar o `lint-imports` no CI.

### 2.5 Cohesion — Coesão das classes

**Como foi usado:**
```bash
cohesion --files Bio/SeqIO/AbiIO.py Bio/SeqIO/AceIO.py --verbose
```

**Como a métrica funciona:** para cada classe, o Cohesion olha os atributos de instância definidos em
`__init__` (via `self.x = ...`) e, para cada método, calcula quantos desses atributos ele referencia. O
"Total" da classe é a média dessa proporção entre os métodos. Quanto mais próximo de 100%, mais os métodos
"trabalham juntos" sobre o mesmo estado (alta coesão); coesão baixa é sinal de classe que faz coisas
desconectadas (candidata a ser dividida).

**Resultado por classe:**

| Classe | Coesão antes | Coesão depois |
|---|---|---|
| `AbiIterator` | 66,67% | 66,67% (sem mudança) |
| `AceIterator` | 50,00% | 50,00% (sem mudança) |

A coesão **não mudou**, e isso é esperado, não é um problema não resolvido: a ferramenta conta o atributo
de classe `modes` (`modes = "b"` / `modes = "t"`) como uma "variável" que nenhum método usa via `self`,
porque `modes` é exigido pela interface `SequenceIterator` do `Bio.SeqIO` e não é específico da instância —
ele já limitava o teto do indicador antes da refatoração e continua limitando depois. Como toda a lógica
extraída virou **funções de módulo** (não novos métodos da classe), a proporção interna de
`AbiIterator`/`AceIterator` — que métodos existem e quais atributos eles tocam — ficou exatamente igual.
Extrair para métodos de instância que não usam `self` teria, na verdade, **piorado** a coesão medida (mais
métodos no denominador sem tocar em atributos), então a decisão de extrair para nível de módulo foi
deliberada para não distorcer essa métrica.

---

## 3. Tabela de refatoração — `Bio/SeqIO/AbiIO.py`

| Classe | Método/trecho | Problema identificado | Métrica antes | Refatoração aplicada | Métrica depois |
|---|---|---|---|---|---|
| `AbiIterator` | `__next__` | Método com complexidade muito alta | 20 | Dividido em 4 funções de módulo, cada uma com uma única responsabilidade: `_read_abi_header` (lê e valida o header), `_extract_abi_tags` (extrai as tags do diretório do arquivo), `_abi_file_stem` (calcula o nome do arquivo sem extensão) e `_build_abi_record` (monta o `SeqRecord` final) — `__next__` passou a só chamar essas 4 funções em sequência | 5 |
| `AbiIterator` | `__next__` | Complexidade cognitiva muito alta | 20 | Mesma divisão em 4 funções descrita acima | 5 |
| `AbiIterator` | classe | Complexidade geral elevada | 13 | Divisão da lógica em funções menores, redistribuindo a responsabilidade entre elas | 6 |
| `AbiIO.py` | arquivo | Baixa manutenibilidade | 49.55 | Extração de métodos e organização das responsabilidades | 52.21 |
| `AbiIterator` | `_parse_tag_data` | Muitos caminhos de execução no método | 12 | Extração da lógica condicional para um dicionário de despacho (`_ELEM_CODE_CONVERTERS`) | 5 |
| `AbiIterator` | `_abi_trim` | Bloco condicional desnecessariamente aninhado | 6 | Reescrita como retorno antecipado (*early return*) | 6 |
| `AbiIterator` | `_abi_parse_header` | Muitas variáveis locais no método | 17 | Remoção de variáveis não utilizadas/desnecessárias | 15 |
| `AbiIO.py` | arquivo | Muitos alertas de code smell | 33 | Extração de métodos e organização das responsabilidades | 1 |

### Métricas analisadas em `AbiIO.py` e onde foram identificadas

| Métrica | Ferramenta | Onde foi identificada | Antes | Depois |
|---|---|---|---|---|
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `AbiIterator.__next__` | 20 | 5 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `_parse_tag_data` | 12 | 5 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `_abi_trim` | 6 | 6 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `AbiIterator` (média da classe) | 13 | 6 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | arquivo inteiro (média entre todas as funções) | 6.86 | 4.09 |
| Complexidade cognitiva | `cognitive-complexity` | `AbiIterator.__next__` | 20 | 5 |
| Complexidade cognitiva | `cognitive-complexity` | `_parse_tag_data` | 16 | 5 |
| Complexidade cognitiva | `cognitive-complexity` | `_abi_trim` | 12 | 6 |
| Complexidade cognitiva | `cognitive-complexity` | arquivo inteiro (média) | 7.86 | 3.09 |
| Índice de manutenibilidade (MI) | Radon (`radon mi`) | arquivo inteiro (métrica não é por método) | 49.55 | 52.21 |
| Code smells | Pylint | arquivo inteiro — a maior parte das ocorrências estava em `__next__` (`too-many-branches`, `too-many-statements`, `too-many-locals`) | 33 | 1 |
| Duplicação de código | Pylint (`R0801`) + leitura manual | bloco `try/except AttributeError` (extração do nome do arquivo) repetido dentro de `__next__` | 1 bloco duplicado | 0 blocos |
| Coesão | Cohesion | classe `AbiIterator` | 66.67% | 66.67% (sem mudança — ver §2.5) |
| Acoplamento | Import Linter | `AbiIterator` — verificação de que não importa `AceIO` diretamente | sem contrato definido | KEPT |

## 4. Tabela de refatoração — `Bio/SeqIO/AceIO.py`

| Classe | Método/trecho | Problema identificado | Métrica antes | Refatoração aplicada | Métrica depois |
|---|---|---|---|---|---|
| `AceIterator` | `__next__` | Método concentrava muitas responsabilidades | 7 | Dividido em 2 funções de módulo, cada uma com uma responsabilidade: `_normalize_consensus_sequence` (normaliza a sequência) e `_consensus_quality_scores` (calcula as qualidades) — `__next__` ficou só chamando as duas | 2 |
| `AceIterator` | `__next__` | Complexidade cognitiva elevada | 7 | Mesma divisão em 2 funções descrita acima | 1 |
| `AceIterator` | classe | Alta complexidade geral da classe | 5 | Divisão das responsabilidades entre funções de módulo | 3 |
| `AceIO.py` | arquivo | Manutenibilidade abaixo do ideal | 67.28 | Extração de métodos e organização das responsabilidades | 67.81 |

### Métricas analisadas em `AceIO.py` e onde foram identificadas

| Métrica | Ferramenta | Onde foi identificada | Antes | Depois |
|---|---|---|---|---|
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `AceIterator.__next__` | 7 | 2 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | `AceIterator` (média da classe) | 5 | 3 |
| Complexidade ciclomática (CC) | Radon (`radon cc`) | arquivo inteiro (média) | 4.00 | 2.50 |
| Complexidade cognitiva | `cognitive-complexity` | `AceIterator.__next__` | 7 | 1 |
| Complexidade cognitiva | `cognitive-complexity` | arquivo inteiro (média) | 3.50 | 1.50 |
| Índice de manutenibilidade (MI) | Radon (`radon mi`) | arquivo inteiro (métrica não é por método) | 67.28 | 67.81 |
| Code smells | Pylint | comentário `# TODO` em `__next__` e nome do módulo `AceIO` — as 2 ocorrências foram mantidas de propósito (ver justificativa abaixo) | 2 | 2 |
| Duplicação de código | Pylint (`R0801`) + leitura manual | nenhum bloco duplicado encontrado | 0 blocos | 0 blocos |
| Coesão | Cohesion | classe `AceIterator` | 50.00% | 50.00% (sem mudança — ver §2.5) |
| Acoplamento | Import Linter | `AceIterator` — verificação de que não importa `AbiIO` diretamente | sem contrato definido | KEPT |

**Itens não alterados de propósito (não são pendências):**
- `_AbiTrimIterator` e o nome do módulo `AceIO` disparam `invalid-name` no Pylint, mas **não foram renomeados** — `Bio/SeqIO/__init__.py` faz o *dispatch* de formato usando esses nomes exatos (`"abi-trim": AbiIO._AbiTrimIterator`, módulo `AceIO` para o formato `"ace"`); renomear quebraria a API pública do Biopython.
- O comentário `# TODO - Supporting reads...` em `AceIterator.__next__` (Pylint `fixme`) foi mantido porque é uma nota de trabalho futuro legítima do autor original do parser, não um code smell.
- A coesão (Cohesion) das duas classes não mudou porque o atributo `modes`, exigido pela interface `SequenceIterator`, não é usado via `self` por nenhum método — isso já limitava o indicador antes da refatoração e continua limitando depois; extrair para métodos de instância (em vez de funções de módulo) teria piorado ainda mais essa métrica. Detalhes em §2.5.

---

## 5. Validação funcional

Não foram criados testes automatizados (conforme pedido), mas antes de fechar a refatoração cada parser
foi comparado, campo a campo, contra uma cópia intocada do Biopython 1.88 (`pip download`), usando os
arquivos de exemplo já existentes em `Tests/Abi/` e `Tests/Ace/` (`310.ab1`, `3730.ab1`, `test.fsa`,
`empty.ab1`, `consed_sample.ace`, `contig1.ace`, `seq.cap.ace`, nos modos `abi`, `abi-trim` e `ace`). Os
resultados (`id`, `seq`, `len`, `run_start`, `phred_quality`) bateram exatamente com a versão original em
todos os casos — a refatoração foi apenas estrutural (extração de funções), sem mudança de comportamento.

## 6. Dashboard (Streamlit + Plotly)

O dashboard existente em `dashboard/dashboard.py` (que já cobria `Bio/Align/bigpsl.py`) foi estendido com
uma segunda seção completa para `AbiIO.py`/`AceIO.py`, reaproveitando as mesmas funções auxiliares
(`criar_grafico_barra`, `executar_comando`, `classe_radon`) e o mesmo estilo visual:

- **Cards** com CC de `__next__` e MI de cada arquivo (antes/depois);
- **Gráfico 1** — complexidade ciclomática de `AbiIterator.__next__` e `AceIterator.__next__` (antes x depois);
- **Gráfico 2** — complexidade ciclomática média do arquivo (antes x depois);
- **Gráfico 3** — Maintainability Index dos dois arquivos lado a lado;
- **Gráfico 4** — coesão das duas classes (antes x depois), com a explicação do artefato do atributo `modes`;
- **Gráfico 5** — complexidade de cada função/método após a refatoração, agrupada por arquivo;
- **Gráfico 6** — status dos dois contratos do Import Linter (`AbiIO → AceIO` e `AceIO → AbiIO`);
- **Gráfico 7** — ocorrências do Pylint antes x depois, com tabela das ocorrências restantes.

Diferente da seção do `bigpsl.py` (que lê valores "depois" de JSONs estáticos), a nova seção executa
**Radon, Pylint, Cohesion e Import Linter ao vivo** (via `subprocess`, com `@st.cache_data`) contra o
código já refatorado em `Bio/SeqIO/AbiIO.py`/`AceIO.py` toda vez que o dashboard roda — só os valores
"antes" ficam fixos como constantes, porque o código anterior à refatoração não existe mais no repositório.
O dashboard foi validado com `streamlit.testing.v1.AppTest` (roda o script sem navegador) e não apresentou
nenhuma exceção.

Para rodar:
```bash
venv/Scripts/streamlit run biopython/dashboard/dashboard.py
```

## 7. Resumo

| Indicador | AbiIO.py — antes | AbiIO.py — depois | AceIO.py — antes | AceIO.py — depois |
|---|---|---|---|---|
| CC de `__next__` | 20 (C) | 5 (A) | 7 (B) | 2 (A) |
| Cognitiva de `__next__` | 20 | 5 | 7 | 1 |
| CC média do arquivo | 6,86 | 4,09 | 4,00 | 2,50 |
| MI do arquivo | 49,55 (A) | 52,21 (A) | 67,28 (A) | 67,81 (A) |
| Ocorrências Pylint | 33 | 1 | 2 | 2 (intencionais) |
| Nota Pylint | 7,92/10 | 9,93/10 | 9,41/10 | 9,46/10 |
| Coesão (Cohesion) | 66,67% | 66,67% | 50,00% | 50,00% |
| Acoplamento (Import Linter) | não medido | KEPT | não medido | KEPT |

---

*Este arquivo não foi commitado — o commit fica a critério do usuário, conforme solicitado.*
