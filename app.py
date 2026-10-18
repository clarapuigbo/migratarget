import streamlit as st
import pandas as pd
import os
import time
import math
import datetime
import re
import requests
from google import genai
from google.genai import types
from Bio import Entrez

st.set_page_config(page_title="MigraTarget", page_icon="🎯", layout="centered")

# ============================================================
# CONFIGURACIÓN
# ============================================================
BASE_URL = "https://api.fda.gov/drug/event.json"
FAERS_API_KEY = st.secrets["FAERS_API_KEY"]
GOOGLE_API_KEY = st.secrets["GOOGLE_API_KEY"]
Entrez.email = st.secrets["NCBI_EMAIL"]
Entrez.api_key = st.secrets["NCBI_API_KEY"]

CARPETA_SALIDA = "resultados"
RUTA_FASE1 = os.path.join(CARPETA_SALIDA, "resultados_migrana_faers.xlsx")
RUTA_FASE2 = os.path.join(CARPETA_SALIDA, "resultados_pubmed_migrana.xlsx")

os.makedirs(CARPETA_SALIDA, exist_ok=True)

CARACTERIZACION_LABELS = {"1": "sospechoso primario", "2": "sospechoso secundario", "3": "concomitante/interactuante"}

# ============================================================
# FASE 1 — FAERS (señales de disproporcionalidad, ROR)
# ============================================================
def pedir_con_reintentos(params, intentos=3):
    for intento in range(1, intentos + 1):
        try:
            return requests.get(BASE_URL, params=params, timeout=60)
        except requests.exceptions.RequestException:
            time.sleep(3)
    return None

def contar(search_query):
    params = {"search": search_query, "limit": 1, "api_key": FAERS_API_KEY}
    r = pedir_con_reintentos(params)
    if r is None or r.status_code != 200:
        return None
    return r.json().get("meta", {}).get("results", {}).get("total", 0)

def calcular_ror(farmaco, a_mas_c, total, reaccion="migraine"):
    a = contar(f'patient.reaction.reactionmeddrapt:"{reaccion}" AND patient.drug.medicinalproduct:"{farmaco}"')
    a_mas_b = contar(f'patient.drug.medicinalproduct:"{farmaco}"')
    if a is None or a_mas_b is None:
        return None
    b = a_mas_b - a
    c = a_mas_c - a
    d = total - a - b - c
    if a == 0 or b == 0 or c == 0 or d == 0:
        return None
    ror = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    ic_bajo = math.exp(math.log(ror) - 1.96 * se)
    ic_alto = math.exp(math.log(ror) + 1.96 * se)
    return {"farmaco": farmaco, "ROR": round(ror, 3), "IC95_bajo": round(ic_bajo, 3), "IC95_alto": round(ic_alto, 3)}

@st.cache_data(ttl=3600)
def obtener_totales_referencia():
    a_mas_c = contar('patient.reaction.reactionmeddrapt:"migraine"')
    total = contar('_exists_:safetyreportid')
    return a_mas_c, total

# Fármacos de mecanismos ya conocidos para migraña/cefalea, según el protocolo
FARMACOS_EXCLUIDOS = [
    "erenumab", "aimovig", "fremanezumab", "ajovy", "galcanezumab", "emgality",
    "eptinezumab", "vyepti", "rimegepant", "nurtec", "ubrogepant", "ubrelvy", "atogepant", "qulipta",
    "sumatriptan", "rizatriptan", "zolmitriptan", "almotriptan", "eletriptan", "frovatriptan", "naratriptan",
    "topiramate", "topamax", "valproate", "divalproex", "depakote", "gabapentin",
    "propranolol", "metoprolol", "timolol", "candesartan", "verapamil", "flunarizine", "amlodipine",
    "botulinum", "botox", "onabotulinumtoxina",
    "ibuprofen", "naproxen", "aspirin", "acetaminophen", "paracetamol", "acetylsalicylic",
    "diclofenac", "ketorolac", "oxycodone", "morphine", "tramadol", "codeine",
    "semaglutide", "ozempic", "wegovy", "liraglutide", "dulaglutide", "tirzepatide",
    "baclofen",
    "magnesium", "riboflavin", "coenzyme q10", "melatonin", "feverfew",
    "vitamin c", "vitamin d", "ascorbic acid", "ascorbate", "cholecalciferol", "ergocalciferol",
    "calcium glubionate", "calcium ascorbate", "folic acid",
]

# Entradas de FAERS que no son fármacos (vacunas, antisépticos) — ruido de la base de datos
NO_FARMACOS = [
    "vaccine", "vaccin", "azd-1222", "comirnaty", "spikevax",
    "hand sanitizer", "sanitizer", "isopropyl alcohol",
]

def es_farmaco_conocido(farmaco):
    nombre = farmaco.lower()
    return any(excluido in nombre for excluido in FARMACOS_EXCLUIDOS)

def es_no_farmacologico(farmaco):
    nombre = farmaco.lower()
    return any(termino in nombre for termino in NO_FARMACOS)

def obtener_clase_atc(farmaco):
    """Devuelve la clase ATC (mecanismo/categoria) de un farmaco via RxNorm+RxClass (NIH, gratis).
    NOTA/limitacion documentada: cuando un farmaco tiene varias clases ATC (distintas formas o usos),
    se coge la primera que devuelve la API -- aproximacion razonable para un grafico agregado,
    no siempre la mas representativa (ej. diclofenac puede salir como "dermatologicals" por su
    forma en gel, en vez de como antiinflamatorio oral)."""
    try:
        r = requests.get(
            "https://rxnav.nlm.nih.gov/REST/rxcui.json",
            params={"name": farmaco, "search": 2}, timeout=10
        )
        rxcui_list = r.json().get("idGroup", {}).get("rxnormId", [])
        if not rxcui_list:
            return "Sin clasificar"
        rxcui = rxcui_list[0]

        r2 = requests.get(
            "https://rxnav.nlm.nih.gov/REST/rxclass/class/byRxcui.json",
            params={"rxcui": rxcui, "relaSource": "ATC"}, timeout=10
        )
        lista = r2.json().get("rxclassDrugInfoList", {}).get("rxclassDrugInfo", [])
        if not lista:
            return "Sin clasificar"
        return lista[0]["rxclassMinConceptItem"]["className"]
    except Exception:
        return "Sin clasificar"

def es_significativo(fila):
    return fila["IC95_bajo"] > 1 or fila["IC95_alto"] < 1

def farmacos_con_caracterizacion(reporte):
    drugs = reporte.get("patient", {}).get("drug", [])
    partes = []
    for d in drugs:
        nombre = d.get("medicinalproduct", "desconocido")
        codigo = d.get("drugcharacterization")
        etiqueta = CARACTERIZACION_LABELS.get(codigo, f"código {codigo}")
        partes.append(f"{nombre} ({etiqueta})")
    return "; ".join(partes)

def actualizar_datos_faers(log=lambda msg: None):
    """Pipeline completo de Fase 1: descarga, filtra, calcula ROR, aplica filtros, guarda Excel."""
    search_query = 'patient.reaction.reactionmeddrapt:"migraine"'
    registros = []
    limit = 100
    skip = 0
    MAX_REGISTROS = 500  # muestra de prueba

    while skip < MAX_REGISTROS:
        log(f"Descargando desde skip={skip}...")
        params = {"search": search_query, "sort": "receivedate:desc", "limit": limit, "skip": skip, "api_key": FAERS_API_KEY}
        r = pedir_con_reintentos(params)
        if r is None or r.status_code != 200:
            log(f"Aviso: fallo al descargar en skip={skip}. Sigo con lo ya descargado.")
            break
        data = r.json()
        if "results" not in data:
            break
        registros.extend(data["results"])
        skip += limit
        if len(data["results"]) < limit:
            break

    def tiene_sospechoso_primario(reporte):
        drugs = reporte.get("patient", {}).get("drug", [])
        return any(d.get("drugcharacterization") == "1" for d in drugs)

    registros_ps = [r for r in registros if tiene_sospechoso_primario(r)]
    df = pd.DataFrame(registros_ps)
    if df.empty:
        return "No se descargaron reportes. Revisa la clave de openFDA."

    df["safetyreportversion"] = pd.to_numeric(df.get("safetyreportversion", 1), errors="coerce").fillna(1)
    df_limpio = (
        df.sort_values("safetyreportversion", ascending=False)
        .drop_duplicates(subset="safetyreportid", keep="first")
        .copy()
    )

    def farmacos_sospechosos_primarios(reporte):
        drugs = reporte.get("patient", {}).get("drug", [])
        return [d.get("medicinalproduct") for d in drugs if d.get("drugcharacterization") == "1"]

    todos_los_farmacos = []
    for _, fila in df_limpio.iterrows():
        todos_los_farmacos.extend(farmacos_sospechosos_primarios(fila))

    tabla_candidatos = pd.Series(todos_los_farmacos).value_counts().reset_index()
    tabla_candidatos.columns = ["farmaco", "num_reportes_migrana_ps"]
    tabla_candidatos = tabla_candidatos[~tabla_candidatos["farmaco"].apply(es_no_farmacologico)]
    candidatos = tabla_candidatos[tabla_candidatos["num_reportes_migrana_ps"] >= 10]

    obtener_totales_referencia.clear()
    a_mas_c, total = obtener_totales_referencia()

    auditoria = []
    for i, farmaco in enumerate(candidatos["farmaco"], start=1):
        log(f"[{i}/{len(candidatos)}] Calculando ROR de {farmaco}...")
        r = calcular_ror(farmaco, a_mas_c, total)

        if r is None:
            auditoria.append({
                "farmaco": farmaco, "ROR": None, "IC95_bajo": None, "IC95_alto": None,
                "incluido_final": False,
                "motivo_descarte": "No se pudo calcular el ROR (datos insuficientes en alguna celda de la tabla 2x2)"
            })
            continue

        significativo = es_significativo(r)
        conocido = es_farmaco_conocido(farmaco)

        if not significativo:
            motivo = "No estadísticamente significativo (el IC95% cruza 1)"
        elif conocido:
            motivo = "Mecanismo farmacológico ya conocido para migraña (excluido por protocolo)"
        else:
            motivo = ""

        auditoria.append({**r, "incluido_final": (significativo and not conocido), "motivo_descarte": motivo})

    if auditoria:
        df_auditoria = pd.DataFrame(auditoria).sort_values("ROR", ascending=False, na_position="last")
    else:
        df_auditoria = pd.DataFrame(columns=["farmaco", "ROR", "IC95_bajo", "IC95_alto", "incluido_final", "motivo_descarte"])

    df_ror = df_auditoria[df_auditoria["ROR"].notna()][["farmaco", "ROR", "IC95_bajo", "IC95_alto"]]
    df_final = df_auditoria[df_auditoria["incluido_final"] == True][["farmaco", "ROR", "IC95_bajo", "IC95_alto"]]

    df_limpio["farmacos_y_rol"] = df_limpio.apply(farmacos_con_caracterizacion, axis=1)
    otras_columnas = [c for c in df_limpio.columns if c != "farmacos_y_rol"]
    df_limpio = df_limpio[["farmacos_y_rol"] + otras_columnas]

    # Clase ATC (RxClass, NIH) solo para los candidatos finales -- son pocos, no todos los evaluados
    df_final = df_final.copy()
    clases_atc = []
    for i, farmaco in enumerate(df_final["farmaco"], start=1):
        log(f"[{i}/{len(df_final)}] Clasificando ATC: {farmaco}...")
        clases_atc.append(obtener_clase_atc(farmaco))
    df_final["clase_atc"] = clases_atc

    with pd.ExcelWriter(RUTA_FASE1) as writer:
        tabla_candidatos.to_excel(writer, sheet_name="candidatos", index=False)
        df_limpio.to_excel(writer, sheet_name="reportes_limpios", index=False)
        df_auditoria.to_excel(writer, sheet_name="auditoria_seleccion", index=False)
        df_ror.to_excel(writer, sheet_name="ROR", index=False)
        df_final.to_excel(writer, sheet_name="ROR_final", index=False)

    return f"Actualizado: {len(df_final)} candidatos finales tras filtros (de {len(df_auditoria)} evaluados)"

# ============================================================
# FASE 2 — PubMed: ¿existe literatura previa? (sin IA, solo código)
# ============================================================
SUFIJOS_QUIMICOS = [
    " hydrochloride", " sodium", " sulfate", " maleate", " acetate", " citrate",
    " besylate", " mesylate", " succinate", " phosphate", " bromide", " tartrate",
    " dihydrate", " anhydrous", " monohydrate", " pivalate", " gluconate",
]

def limpiar_nombre_farmaco(farmaco):
    """Quita sales químicas y paréntesis, para que la búsqueda en PubMed encuentre
    resultados aunque el nombre de FAERS venga con la forma química completa."""
    nombre = farmaco.lower()
    for sufijo in SUFIJOS_QUIMICOS:
        if nombre.endswith(sufijo):
            nombre = nombre[: -len(sufijo)]
    nombre = nombre.replace("(", "").replace(")", "")
    return nombre.strip()

def _es_libro(articulo):
    """PubMed indexa 'libros' (StatPearls, capitulos de referencia) con una
    estructura distinta ('BookDocument') a la de los articulos normales
    ('MedlineCitation'). Hay que distinguirlos para leer cada uno bien."""
    return "BookDocument" in articulo

def extraer_pmid(articulo):
    try:
        contenedor = "BookDocument" if _es_libro(articulo) else "MedlineCitation"
        return str(articulo[contenedor]["PMID"])
    except (KeyError, TypeError):
        return "desconocido"

def extraer_titulo(articulo):
    try:
        if _es_libro(articulo):
            return str(articulo["BookDocument"].get("ArticleTitle", "(sin titulo)"))
        return str(articulo["MedlineCitation"]["Article"]["ArticleTitle"])
    except (KeyError, TypeError):
        return "(sin titulo)"

def extraer_fecha_publicacion(articulo):
    """Devuelve la fecha de publicacion como texto (año, o año-mes-dia si esta disponible).
    Funciona tanto para articulos normales como para 'libros' (StatPearls, etc.)."""
    try:
        if _es_libro(articulo):
            pub_date = articulo["BookDocument"]["Book"].get("PubDate", {})
        else:
            pub_date = articulo["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
        anio = pub_date.get("Year", "")
        mes = pub_date.get("Month", "")
        dia = pub_date.get("Day", "")
        if anio:
            return "-".join(p for p in [anio, mes, dia] if p)
        return str(pub_date.get("MedlineDate", "fecha desconocida"))
    except (KeyError, TypeError):
        return "fecha desconocida"

def extraer_autores(articulo):
    """Devuelve 'Apellido Iniciales, Apellido Iniciales, ...' de la lista de autores.
    Funciona tanto para articulos normales como para 'libros' (StatPearls, etc.)."""
    try:
        if _es_libro(articulo):
            lista = articulo["BookDocument"].get("AuthorList", [])
            if not lista:
                lista = articulo["BookDocument"]["Book"].get("AuthorList", [])
        else:
            lista = articulo["MedlineCitation"]["Article"].get("AuthorList", [])
        nombres = []
        for a in lista:
            apellido = a.get("LastName", "")
            inicial = a.get("Initials", "")
            if apellido:
                nombres.append(f"{apellido} {inicial}".strip())
        return ", ".join(nombres) if nombres else "autores no disponibles"
    except (KeyError, TypeError):
        return "autores no disponibles"

def anio_de_fecha(fecha_texto):
    """Extrae el año (4 digitos) de un texto de fecha, para poder ordenar. 0 si no se encuentra."""
    m = re.search(r"\d{4}", str(fecha_texto))
    return int(m.group()) if m else 0

def buscar_evidencia_pubmed(farmaco, retmax=10, intentos=3):
    """Busca si existe literatura previa en PubMed sobre un fármaco y la migraña.
    Devuelve título, fecha, autores y enlace directo de cada artículo encontrado."""
    nombre_busqueda = limpiar_nombre_farmaco(farmaco)
    query = f'("{nombre_busqueda}") AND (Migraine*)'

    for intento in range(1, intentos + 1):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort="pub_date")
            resultado = Entrez.read(handle)
            handle.close()
            total = int(resultado["Count"])
            ids = resultado["IdList"]
            if total == 0:
                return {"estado": "No DATA available", "num_articulos": 0, "articulos": []}

            handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
            articulos_raw = Entrez.read(handle)
            handle.close()

            # Se combinan articulos normales y 'libros' (StatPearls, etc.) -- ambos
            # pueden venir en la respuesta, con estructuras distintas.
            todos_los_articulos_raw = list(articulos_raw.get("PubmedArticle", [])) + list(articulos_raw.get("PubmedBookArticle", []))

            articulos = []
            for articulo in todos_los_articulos_raw:
                pmid = extraer_pmid(articulo)
                articulos.append({
                    "titulo": extraer_titulo(articulo),
                    "fecha_publicacion": extraer_fecha_publicacion(articulo),
                    "autores": extraer_autores(articulo),
                    "pmid": pmid,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                })
            return {"estado": "con evidencia", "num_articulos": total, "articulos": articulos}
        except Exception:
            time.sleep(3 * intento)

    return {"estado": "Error de conexión", "num_articulos": 0, "articulos": []}

def actualizar_datos_pubmed(log=lambda msg: None):
    """Para cada candidato final de Fase 1, comprueba si hay literatura previa en PubMed."""
    if not os.path.exists(RUTA_FASE1):
        return "Error: no hay datos de Fase 1 todavía. Actualiza FAERS primero."

    candidatos = pd.read_excel(RUTA_FASE1, sheet_name="ROR_final")
    if candidatos.empty:
        return "No hay candidatos en Fase 1 para buscar evidencia."

    evidencia_pubmed = []
    for i, farmaco in enumerate(candidatos["farmaco"], start=1):
        log(f"[{i}/{len(candidatos)}] Buscando en PubMed: {farmaco}...")
        resultado = buscar_evidencia_pubmed(farmaco)
        evidencia_pubmed.append({"farmaco": farmaco, **resultado})
        time.sleep(0.4)

    # Ordenado por el articulo mas reciente de cada farmaco (descendente)
    filas_resumen = []
    filas_articulos = []
    for e in evidencia_pubmed:
        anios = [anio_de_fecha(a["fecha_publicacion"]) for a in e["articulos"]]
        anio_mas_reciente = max(anios) if anios else 0
        filas_resumen.append({
            "farmaco": e["farmaco"],
            "num_articulos": e["num_articulos"],
            "estado": e["estado"],
            "anio_mas_reciente": anio_mas_reciente if anio_mas_reciente else "N/A",
        })
        for articulo in e["articulos"]:
            filas_articulos.append({
                "farmaco": e["farmaco"],
                "titulo": articulo["titulo"],
                "fecha_publicacion": articulo["fecha_publicacion"],
                "autores": articulo["autores"],
                "url": articulo["url"],
                "_anio_orden": anio_de_fecha(articulo["fecha_publicacion"]),
            })

    tabla_resumen = pd.DataFrame(filas_resumen).sort_values(
        "anio_mas_reciente", ascending=False, key=lambda s: pd.to_numeric(s, errors="coerce").fillna(0)
    )
    tabla_articulos_encontrados = (
        pd.DataFrame(filas_articulos)
        .sort_values("_anio_orden", ascending=False)
        .drop(columns="_anio_orden")
    )

    with pd.ExcelWriter(RUTA_FASE2) as writer:
        tabla_resumen.to_excel(writer, sheet_name="resumen", index=False)
        tabla_articulos_encontrados.to_excel(writer, sheet_name="articulos_encontrados", index=False)

    return f"Actualizado: {len(candidatos)} candidatos revisados en PubMed"

# ============================================================
# HERRAMIENTAS DEL AGENTE DE CHAT
# ============================================================
a_mas_c_total, total_base = obtener_totales_referencia()

def herramienta_calcular_ror(farmaco: str) -> dict:
    """Calcula el ROR (Reporting Odds Ratio) de un farmaco respecto a reportes de migrana en FAERS.

    Args:
        farmaco: nombre generico del farmaco en ingles, tal como aparece en FAERS (ej. "methotrexate").
    """
    resultado = calcular_ror(farmaco, a_mas_c_total, total_base)
    if resultado is None:
        return {"error": f"No se pudo calcular el ROR para {farmaco} (datos insuficientes)"}
    return resultado

def herramienta_buscar_pubmed(farmaco: str) -> dict:
    """Busca si hay evidencia cientifica en PubMed que relacione un farmaco con la migrana,
    e incluye el enlace directo a cada articulo encontrado para que el usuario pueda abrirlo.

    Args:
        farmaco: nombre generico del farmaco en ingles (ej. "infliximab").
    """
    resultado = buscar_evidencia_pubmed(farmaco)
    if resultado["estado"] != "con evidencia":
        return {"estado": resultado["estado"]}
    articulos_resumidos = [{"titulo": a["titulo"], "url": a["url"]} for a in resultado["articulos"]]
    return {"estado": "con evidencia", "num_articulos": resultado["num_articulos"], "articulos": articulos_resumidos}

@st.cache_resource
def obtener_cliente():
    return genai.Client(api_key=GOOGLE_API_KEY)

client = obtener_cliente()

def preguntar_agente(pregunta):
    respuesta = client.models.generate_content(
        model="gemini-3.6-flash",
        contents=pregunta,
        config=types.GenerateContentConfig(tools=[herramienta_calcular_ror, herramienta_buscar_pubmed])
    )
    return respuesta.text

def fecha_ultima_modificacion(ruta):
    if not os.path.exists(ruta):
        return None
    return datetime.datetime.fromtimestamp(os.path.getmtime(ruta)).strftime("%d/%m/%Y %H:%M")

def mostrar_tabla_bonita(df, alto="420px"):
    """Muestra una tabla con el header en negrita/color de marca -- st.dataframe no lo permite
    (su header se dibuja en un lienzo, no como texto normal), así que se renderiza como HTML."""
    html = df.to_html(classes="tabla-bonita", escape=False, index=False)
    st.markdown(
        f'''<style>
.tabla-bonita {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
.tabla-bonita th {{ background-color: #0A2A2B; color: #FF7A5C; font-weight: 700;
                     padding: 8px; text-align: left; position: sticky; top: 0; }}
.tabla-bonita td {{ padding: 6px 8px; border-bottom: 1px solid #DCF0EE; }}
</style>
<div style="max-height:{alto}; overflow:auto;">{html}</div>''',
        unsafe_allow_html=True,
    )

# ============================================================
# INTERFAZ
# ============================================================
col_logo, col_titulo = st.columns([1, 5])
with col_logo:
    st.image("logo_migratarget.svg", width=90)
with col_titulo:
    st.title("MigraTarget")
    st.caption("Señales FAERS y literatura previa sobre migraña")

st.subheader("Pregunta al agente")
pregunta = st.text_input("¿Qué quieres saber? (ROR de un fármaco, o si hay literatura en PubMed)")
if pregunta:
    with st.spinner("Pensando..."):
        respuesta = preguntar_agente(pregunta)
    st.write(respuesta)

st.divider()
st.subheader("Mantenimiento de datos — Fase 1 (FAERS)")
if st.button("Actualizar candidatos desde FAERS (tarda varios minutos)"):
    estado = st.empty()
    resultado = actualizar_datos_faers(log=lambda msg: estado.text(msg))
    st.success(resultado)

st.subheader("Mantenimiento de datos — Fase 2 (comprobar literatura en PubMed)")
st.caption("Usa los candidatos ya guardados en Fase 1. Solo comprueba si hay artículos o no — sin evaluación por IA.")
if st.button("Comprobar literatura en PubMed"):
    estado = st.empty()
    resultado = actualizar_datos_pubmed(log=lambda msg: estado.text(msg))
    st.success(resultado)

st.divider()
st.subheader("Resultados guardados — Fase 1 (FAERS)")
fecha1 = fecha_ultima_modificacion(RUTA_FASE1)
if fecha1:
    st.caption(f"Última actualización: {fecha1}")
    df_candidatos = pd.read_excel(RUTA_FASE1, sheet_name="ROR_final")
    mostrar_tabla_bonita(df_candidatos)

    if "clase_atc" in df_candidatos.columns:
        st.caption("Distribución por clase ATC (RxClass, NIH) -- aproximación, ver nota en el código")
        conteo_atc = df_candidatos["clase_atc"].value_counts()
        st.bar_chart(conteo_atc)

    with open(RUTA_FASE1, "rb") as f:
        st.download_button("Descargar Excel de Fase 1", data=f, file_name="resultados_migrana_faers.xlsx")
else:
    st.info("Todavía no se ha generado el Excel de Fase 1.")

st.divider()
st.subheader("Trazabilidad de la selección — Fase 1")
if fecha1:
    df_auditoria = pd.read_excel(RUTA_FASE1, sheet_name="auditoria_seleccion")
    df_auditoria["motivo_descarte"] = df_auditoria["motivo_descarte"].fillna("")
    incluidos = int(df_auditoria["incluido_final"].sum())
    st.caption(f"{incluidos} de {len(df_auditoria)} candidatos pasaron todos los filtros")
    resumen_motivos = df_auditoria[df_auditoria["motivo_descarte"] != ""]["motivo_descarte"].value_counts()
    if not resumen_motivos.empty:
        st.caption("Motivos de descarte:")
        mostrar_tabla_bonita(resumen_motivos.reset_index().rename(columns={"index": "motivo", "motivo_descarte": "num_farmacos"}), alto="200px")
    df_auditoria_mostrar = df_auditoria.copy()
    df_auditoria_mostrar["incluido_final"] = df_auditoria_mostrar["incluido_final"].apply(
        lambda v: '<span style="color:#FF7A5C;font-weight:700">Sí</span>' if v
        else '<span style="color:#0A2A2B;font-weight:700">No</span>'
    )
    mostrar_tabla_bonita(df_auditoria_mostrar)
else:
    st.info("Todavía no hay tabla de auditoría — actualiza Fase 1 primero.")

st.divider()
st.subheader("Resultados guardados — Fase 2 (literatura en PubMed)")
fecha2 = fecha_ultima_modificacion(RUTA_FASE2)
if fecha2:
    st.caption(f"Última actualización: {fecha2}")
    st.caption("Ordenado por el artículo más reciente de cada fármaco (año más nuevo primero).")
    df_resumen2 = pd.read_excel(RUTA_FASE2, sheet_name="resumen")
    mostrar_tabla_bonita(df_resumen2)

    with st.expander("Ver artículos encontrados, con enlaces"):
        df_articulos2 = pd.read_excel(RUTA_FASE2, sheet_name="articulos_encontrados")
        df_articulos2["url"] = df_articulos2["url"].apply(lambda u: f'<a href="{u}" target="_blank">Abrir</a>')
        mostrar_tabla_bonita(df_articulos2)

    with open(RUTA_FASE2, "rb") as f:
        st.download_button("Descargar Excel de Fase 2", data=f, file_name="resultados_pubmed_migrana.xlsx")
else:
    st.info("Todavía no se ha generado el Excel de Fase 2.")
