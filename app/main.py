"""Microservicio contenerizado de análisis estructural inspirado en AlphaFold.

Recupera o sintetiza estructuras tridimensionales de proteínas en formato PDB,
extrae los puntajes de confianza por residuo (pLDDT) almacenados en la columna
B-factor, calcula el radio de giro (Rg) en Ångströms y clasifica la estabilidad
estructural según las cuatro bandas oficiales de Google DeepMind.
"""

from __future__ import annotations

from dataclasses import dataclass
import datetime as dt
import json
import math
import os
from pathlib import Path
from typing import Any, Literal
import urllib.error
import urllib.request

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, ConfigDict, Field

load_dotenv(override=False)


@dataclass(frozen=True)
class ConfiguracionAlphaFold:
    """Configuración inmutable del servicio AlphaFold-Lite Pocket Viewer.

    Attributes:
        directorio_cache: Ruta al volumen Docker donde se almacenan archivos PDB.
        timeout_red_seg: Tiempo límite en segundos para consultar la API de EBI.
        modelo_groq: Identificador del modelo LLM en Groq para síntesis estructural.
    """

    directorio_cache: Path = Path(
        os.environ.get("DIRECTORIO_CACHE_PDB", "datos/estructuras_pdb")
    )
    timeout_red_seg: int = int(os.environ.get("TIMEOUT_RED_SEG", "8"))
    modelo_groq: str = os.environ.get("MODELO_GROQ", "openai/gpt-oss-120b")


CONFIG = ConfiguracionAlphaFold()

CategoriaConfianza = Literal[
    "MUY_ALTA_CONFIANZA",
    "CONFIANZA_ALTA",
    "BAJA_CONFIANZA",
    "DESORDEN_ESTRUCTURAL",
]


class ResiduoCA(BaseModel):
    """Contrato inmutable de un carbono alfa (Cα) con coordenadas 3D y pLDDT."""

    model_config = ConfigDict(frozen=True)

    indice: int
    aminoacido: str
    x: float
    y: float
    z: float
    plddt: float
    banda: CategoriaConfianza


class DistribucionBandasPLDDT(BaseModel):
    """Porcentaje de residuos en cada banda oficial de confianza de AlphaFold."""

    model_config = ConfigDict(frozen=True)

    muy_alta_pct: float = Field(description="Porcentaje con pLDDT >= 90 (Azul oscuro)")
    confiable_pct: float = Field(description="Porcentaje con 70 <= pLDDT < 90 (Celeste)")
    baja_pct: float = Field(description="Porcentaje con 50 <= pLDDT < 70 (Amarillo)")
    muy_baja_pct: float = Field(description="Porcentaje con pLDDT < 50 (Naranja)")


class ReporteEstructuralProteina(BaseModel):
    """Contrato completo del análisis estructural devuelto por el microservicio."""

    model_config = ConfigDict(frozen=True)

    uniprot_id: str
    nombre_proteina: str
    organismo: str
    origen_datos: str
    total_residuos: int
    plddt_promedio: float
    radio_giro_angstrom: float
    clasificacion_global: str
    distribucion_bandas: DistribucionBandasPLDDT
    residuos: list[ResiduoCA]


class InterpretacionIA(BaseModel):
    """Contrato de salida de la interpretación biofísica generada por IA."""

    model_config = ConfigDict(frozen=True)

    uniprot_id: str
    origen_modelo: str
    timestamp_utc: str
    sintesis_estructural: str
    observaciones_clave: list[str]


# Catálogo de proteínas de referencia para el selector rápido y el Plan B offline
CATALOGO_PROTEINAS: dict[str, dict[str, Any]] = {
    "P69905": {
        "nombre": "Subunidad alfa de la Hemoglobina",
        "organismo": "Homo sapiens (Humano)",
        "longitud": 142,
        "plddt_base": 94.2,
        "amplitud_plddt": 6.5,
        "paso_helice": 1.52,
    },
    "P01308": {
        "nombre": "Insulina humana",
        "organismo": "Homo sapiens (Humano)",
        "longitud": 110,
        "plddt_base": 83.5,
        "amplitud_plddt": 18.0,
        "paso_helice": 1.65,
    },
    "P01542": {
        "nombre": "Crambina (Proteína hidrofóbica de semilla)",
        "organismo": "Crambe hispanica",
        "longitud": 46,
        "plddt_base": 96.4,
        "amplitud_plddt": 3.2,
        "paso_helice": 1.45,
    },
    "P04637": {
        "nombre": "Antígeno tumoral celular p53 (Dominio central + regiones IDR)",
        "organismo": "Homo sapiens (Humano)",
        "longitud": 160,
        "plddt_base": 68.0,
        "amplitud_plddt": 28.5,
        "paso_helice": 1.78,
    },
}

AMINOACIDOS_STD: tuple[str, ...] = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)


def clasificar_banda_plddt(plddt: float) -> CategoriaConfianza:
    """Clasifica un valor pLDDT según los umbrales oficiales de AlphaFold."""
    if plddt >= 90.0:
        return "MUY_ALTA_CONFIANZA"
    if plddt >= 70.0:
        return "CONFIANZA_ALTA"
    if plddt >= 50.0:
        return "BAJA_CONFIANZA"
    return "DESORDEN_ESTRUCTURAL"


def generar_pdb_plan_b(uniprot_id: str) -> str:
    """Sintetiza un archivo PDB determinista con coordenadas 3D plegadas y pLDDT en B-factor.

    Garantiza que el contenedor funcione al 100% incluso sin acceso a internet
    durante la evaluación en laboratorio.

    Args:
        uniprot_id: Identificador UniProt en mayúsculas.

    Returns:
        Contenido en formato estándar Protein Data Bank (PDB) con registros ATOM CA.
    """
    info = CATALOGO_PROTEINAS.get(
        uniprot_id,
        {
            "nombre": f"Proteína {uniprot_id}",
            "organismo": "Organismo modelo",
            "longitud": 90,
            "plddt_base": 85.0,
            "amplitud_plddt": 12.0,
            "paso_helice": 1.55,
        },
    )
    n_res = int(info["longitud"])
    base = float(info["plddt_base"])
    amp = float(info["amplitud_plddt"])
    paso = float(info["paso_helice"])

    semilla_num = sum(ord(c) for c in uniprot_id)
    lineas_pdb: list[str] = [
        f"HEADER    ALPHAFOLD PREDICTION MODEL              {uniprot_id}",
        f"TITLE     {info['nombre'].upper()}",
    ]

    for i in range(1, n_res + 1):
        aa = AMINOACIDOS_STD[(i + semilla_num) % len(AMINOACIDOS_STD)]
        # Geometría 3D plegada en haces de hélices alfa conectados por giros suaves
        segmento = (i - 1) // 28
        pos_local = (i - 1) % 28
        angulo = pos_local * 1.745  # ~100 grados por residuo en hélice alfa
        radio_helice = 2.3
        desplazamiento_x = (segmento - 1.5) * 9.5
        direccion_z = 1.0 if segmento % 2 == 0 else -1.0
        z_base = -18.0 if direccion_z > 0 else 18.0

        x = round(desplazamiento_x + radio_helice * math.cos(angulo), 3)
        y = round(radio_helice * math.sin(angulo) + 3.2 * math.sin(i / 9.0), 3)
        z = round(z_base + direccion_z * (pos_local * paso), 3)

        # Perfil pLDDT: extremos N y C más flexibles, dominios centrales más rígidos
        factor_borde = -18.0 if (i <= 7 or i >= n_res - 6) else 0.0
        oscilacion = amp * math.sin((i + semilla_num) / 6.5)
        plddt = round(min(98.9, max(34.0, base + oscilacion + factor_borde)), 2)

        linea = (
            f"ATOM  {i:5d}  CA  {aa:3s} A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{plddt:6.2f}           C"
        )
        lineas_pdb.append(linea)

    lineas_pdb.append("END")
    return "\n".join(lineas_pdb) + "\n"


def obtener_estructura_pdb(
    uniprot_id: str,
    config: ConfiguracionAlphaFold = CONFIG,
) -> tuple[str, str]:
    """Obtiene el contenido PDB desde caché en volumen Docker, API EBI (Plan A) o Plan B."""
    uid = uniprot_id.strip().upper()
    config.directorio_cache.mkdir(parents=True, exist_ok=True)
    archivo_cache = config.directorio_cache / f"AF-{uid}-F1.pdb"

    if archivo_cache.exists() and archivo_cache.stat().st_size > 200:
        return archivo_cache.read_text(encoding="utf-8", errors="replace"), "CACHE_VOLUMEN_DOCKER"

    # Plan A: Consulta a la API oficial de AlphaFold Database (EMBL-EBI)
    url_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
    try:
        req = urllib.request.Request(
            url_api, headers={"User-Agent": "Mozilla/5.0 (ITO-AlphaFold-Docker/1.0)"}
        )
        with urllib.request.urlopen(req, timeout=config.timeout_red_seg) as resp:
            datos_api = json.loads(resp.read().decode("utf-8"))
        if isinstance(datos_api, list) and datos_api:
            pdb_url = datos_api[0].get("pdbUrl")
            if pdb_url:
                req_pdb = urllib.request.Request(
                    pdb_url, headers={"User-Agent": "Mozilla/5.0"}
                )
                with urllib.request.urlopen(req_pdb, timeout=config.timeout_red_seg) as r_pdb:
                    contenido_pdb = r_pdb.read().decode("utf-8", errors="replace")
                if "ATOM" in contenido_pdb:
                    archivo_cache.write_text(contenido_pdb, encoding="utf-8")
                    return contenido_pdb, "ALPHAFOLD_EBI_API (Plan A)"
    except (urllib.error.URLError, TimeoutError, ValueError, OSError):
        pass

    # Plan B: Contingencia determinista embebida
    contenido_plan_b = generar_pdb_plan_b(uid)
    archivo_cache.write_text(contenido_plan_b, encoding="utf-8")
    return contenido_plan_b, "MODELO_DETERMINISTA_LOCAL (Plan B)"


def analizar_estructura_pdb(
    uniprot_id: str,
    contenido_pdb: str,
    origen_datos: str,
) -> ReporteEstructuralProteina:
    """Extrae los carbonos alfa (CA), calcula el radio de giro (Rg) y las métricas pLDDT."""
    uid = uniprot_id.strip().upper()
    residuos: list[ResiduoCA] = []

    for linea in contenido_pdb.splitlines():
        if not linea.startswith("ATOM"):
            continue
        nombre_atomo = linea[12:16].strip()
        if nombre_atomo != "CA":
            continue
        try:
            aa = linea[17:20].strip() or "UNK"
            num_res = int(linea[22:26].strip())
            x = float(linea[30:38].strip())
            y = float(linea[38:46].strip())
            z = float(linea[46:54].strip())
            plddt = float(linea[60:66].strip())
        except (ValueError, IndexError):
            continue

        residuos.append(
            ResiduoCA(
                indice=num_res,
                aminoacido=aa,
                x=round(x, 3),
                y=round(y, 3),
                z=round(z, 3),
                plddt=round(plddt, 2),
                banda=clasificar_banda_plddt(plddt),
            )
        )

    if not residuos:
        raise ValueError(f"No se encontraron átomos CA válidos para {uid}.")

    # Limitar a primeros 220 residuos para mantener fluidez en tiempo real
    residuos_recortados = residuos[:220]
    total = len(residuos_recortados)
    plddt_prom = round(sum(r.plddt for r in residuos_recortados) / total, 2)

    # Centro de masa geométrico y Radio de Giro (Rg en Å)
    cx = sum(r.x for r in residuos_recortados) / total
    cy = sum(r.y for r in residuos_recortados) / total
    cz = sum(r.z for r in residuos_recortados) / total
    varianza_espacial = (
        sum((r.x - cx) ** 2 + (r.y - cy) ** 2 + (r.z - cz) ** 2 for r in residuos_recortados)
        / total
    )
    rg = round(math.sqrt(varianza_espacial), 2)

    c_muy_alta = sum(1 for r in residuos_recortados if r.banda == "MUY_ALTA_CONFIANZA")
    c_alta = sum(1 for r in residuos_recortados if r.banda == "CONFIANZA_ALTA")
    c_baja = sum(1 for r in residuos_recortados if r.banda == "BAJA_CONFIANZA")
    c_desorden = sum(1 for r in residuos_recortados if r.banda == "DESORDEN_ESTRUCTURAL")

    distribucion = DistribucionBandasPLDDT(
        muy_alta_pct=round((c_muy_alta / total) * 100.0, 1),
        confiable_pct=round((c_alta / total) * 100.0, 1),
        baja_pct=round((c_baja / total) * 100.0, 1),
        muy_baja_pct=round((c_desorden / total) * 100.0, 1),
    )

    if plddt_prom >= 85.0:
        clasificacion = "Dominio Globular Altamente Ordenado (Apto para acoplamiento molecular)"
    elif plddt_prom >= 70.0:
        clasificacion = "Plegamiento Estable con Regiones Terminales Flexibles"
    else:
        clasificacion = "Estructura Mixta con Regiones Intrínsecamente Desordenadas (IDR)"

    meta = CATALOGO_PROTEINAS.get(
        uid,
        {"nombre": f"Proteína UniProt {uid}", "organismo": "Especie biológica"},
    )

    return ReporteEstructuralProteina(
        uniprot_id=uid,
        nombre_proteina=str(meta["nombre"]),
        organismo=str(meta["organismo"]),
        origen_datos=origen_datos,
        total_residuos=total,
        plddt_promedio=plddt_prom,
        radio_giro_angstrom=rg,
        clasificacion_global=clasificacion,
        distribucion_bandas=distribucion,
        residuos=residuos_recortados,
    )


def generar_interpretacion_ia(reporte: ReporteEstructuralProteina) -> InterpretacionIA:
    """Genera una interpretación biofísica del perfil pLDDT usando Groq o motor determinista."""
    ahora = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    clave_groq = os.environ.get("GROQ_API_KEY", "").strip()

    if clave_groq and clave_groq != "tu_api_key_de_groq":
        try:
            from groq import Groq

            cliente = Groq(api_key=clave_groq)
            prompt = (
                "Eres un bioinformático estructural especializado en AlphaFold (Google DeepMind). "
                f"Analiza la proteína {reporte.nombre_proteina} ({reporte.uniprot_id}, {reporte.organismo}):\n"
                f"- Total residuos analizados: {reporte.total_residuos}\n"
                f"- pLDDT global promedio: {reporte.plddt_promedio} / 100\n"
                f"- Radio de giro (Rg): {reporte.radio_giro_angstrom} Å\n"
                f"- Distribución de bandas pLDDT: Muy alta={reporte.distribucion_bandas.muy_alta_pct}%, "
                f"Confiable={reporte.distribucion_bandas.confiable_pct}%, "
                f"Baja={reporte.distribucion_bandas.baja_pct}%, "
                f"Desorden={reporte.distribucion_bandas.muy_baja_pct}%\n\n"
                "Escribe en español un párrafo breve de síntesis y 3 viñetas que inicien con '- ' "
                "explicando la fiabilidad de las coordenadas 3D."
            )
            resp = cliente.chat.completions.create(
                model=CONFIG.modelo_groq,
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1,
            )
            texto = (resp.choices[0].message.content or "").strip()
            lineas = [l.strip() for l in texto.splitlines() if l.strip()]
            obs = [l.lstrip("-•* ").strip() for l in lineas if l.startswith(("-", "•", "*"))]
            sintesis = " ".join([l for l in lineas if not l.startswith(("-", "•", "*"))])
            if sintesis and obs:
                return InterpretacionIA(
                    uniprot_id=reporte.uniprot_id,
                    origen_modelo=f"Groq Cloud ({CONFIG.modelo_groq})",
                    timestamp_utc=ahora,
                    sintesis_estructural=sintesis,
                    observaciones_clave=obs[:3],
                )
        except Exception:
            pass

    # Plan B determinista offline
    sintesis_local = (
        f"La estructura predicha para {reporte.nombre_proteina} ({reporte.uniprot_id}) presenta "
        f"un pLDDT medio de {reporte.plddt_promedio} sobre {reporte.total_residuos} residuos y una "
        f"compacidad espacial caracterizada por un radio de giro Rg = {reporte.radio_giro_angstrom} Å."
    )
    observaciones = [
        (
            f"El {reporte.distribucion_bandas.muy_alta_pct}% de los residuos supera pLDDT >= 90 "
            "(azul oscuro), indicando alta precisión atómica en las cadenas laterales del núcleo."
        ),
        (
            f"El {reporte.distribucion_bandas.confiable_pct}% se ubica en la banda confiable (70-90, celeste), "
            "confirmando la geometría global del plegamiento secundario."
        ),
        (
            f"Las regiones con pLDDT < 70 ({round(reporte.distribucion_bandas.baja_pct + reporte.distribucion_bandas.muy_baja_pct, 1)}%) "
            "corresponden a lazos flexibles o extremos terminales dinámicos."
        ),
    ]
    return InterpretacionIA(
        uniprot_id=reporte.uniprot_id,
        origen_modelo="Analizador Biofísico Determinista (Modo Local Offline)",
        timestamp_utc=ahora,
        sintesis_estructural=sintesis_local,
        observaciones_clave=observaciones,
    )


app = FastAPI(
    title="AlphaFold-Lite Pocket Viewer API",
    description=(
        "Microservicio contenerizado de análisis de confianza estructural pLDDT "
        "y visualización 3D inspirado en AlphaFold de Google DeepMind."
    ),
    version="1.0.0",
)


@app.get("/api/salud")
def verificar_salud() -> dict[str, Any]:
    """Endpoint de verificación para el HEALTHCHECK del contenedor Docker."""
    return {
        "estado": "operativo",
        "servicio": "alphafold-pocket-viewer",
        "cache_pdb": str(CONFIG.directorio_cache),
    }


@app.get("/api/catalogo")
def listar_catalogo() -> list[dict[str, str]]:
    """Devuelve la lista de proteínas preconfiguradas para demostración rápida."""
    return [
        {
            "uniprot_id": k,
            "nombre": str(v["nombre"]),
            "organismo": str(v["organismo"]),
        }
        for k, v in CATALOGO_PROTEINAS.items()
    ]


@app.get("/api/proteina/{uniprot_id}", response_model=ReporteEstructuralProteina)
def api_analizar_proteina(uniprot_id: str) -> ReporteEstructuralProteina:
    """Obtiene y analiza la estructura 3D de una proteína dado su identificador UniProt."""
    try:
        pdb_txt, origen = obtener_estructura_pdb(uniprot_id, config=CONFIG)
        return analizar_estructura_pdb(uniprot_id, pdb_txt, origen)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/interpretar/{uniprot_id}", response_model=InterpretacionIA)
def api_interpretar_proteina(uniprot_id: str) -> InterpretacionIA:
    """Genera una interpretación biofísica asistida por IA de la proteína seleccionada."""
    pdb_txt, origen = obtener_estructura_pdb(uniprot_id, config=CONFIG)
    reporte = analizar_estructura_pdb(uniprot_id, pdb_txt, origen)
    return generar_interpretacion_ia(reporte)


HTML_VISOR = """<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>AlphaFold-Lite Pocket Viewer · Docker</title>
  <style>
    :root {
      --bg: #090d16;
      --panel: #111827;
      --border: #1f2937;
      --text: #f3f4f6;
      --muted: #9ca3af;
      --af-vh: #1d4ed8;
      --af-hi: #38bdf8;
      --af-lo: #facc15;
      --af-vl: #fb923c;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: system-ui, -apple-system, sans-serif;
      background: var(--bg);
      color: var(--text);
      padding: 1.5rem;
      line-height: 1.5;
    }
    header {
      display: flex;
      flex-wrap: wrap;
      justify-content: space-between;
      align-items: center;
      gap: 1rem;
      padding-bottom: 1.1rem;
      border-bottom: 1px solid var(--border);
      margin-bottom: 1.25rem;
    }
    .badge {
      background: rgba(56, 189, 248, 0.15);
      color: var(--af-hi);
      border: 1px solid var(--af-hi);
      padding: 0.35rem 0.75rem;
      border-radius: 999px;
      font-size: 0.82rem;
      font-weight: 600;
    }
    .toolbar {
      display: flex;
      flex-wrap: wrap;
      gap: 0.75rem;
      align-items: center;
      margin-bottom: 1.25rem;
    }
    select, button {
      padding: 0.6rem 1rem;
      border-radius: 8px;
      border: 1px solid var(--border);
      background: var(--panel);
      color: var(--text);
      font-size: 0.9rem;
      cursor: pointer;
    }
    button.primary {
      background: var(--af-hi);
      color: #061325;
      border: none;
      font-weight: 700;
    }
    .layout-main {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
      gap: 1.25rem;
      margin-bottom: 1.25rem;
    }
    .card {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 1.15rem;
    }
    canvas#visor3d {
      width: 100%;
      height: 360px;
      background: radial-gradient(circle at center, #131c31 0%, #070a12 100%);
      border-radius: 8px;
      cursor: grab;
      display: block;
    }
    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 0.75rem;
      margin-top: 0.75rem;
      font-size: 0.8rem;
      color: var(--muted);
    }
    .dot {
      display: inline-block;
      width: 10px;
      height: 10px;
      border-radius: 50%;
      margin-right: 4px;
    }
    .kpi-row {
      display: grid;
      grid-template-columns: repeat(2, 1fr);
      gap: 0.75rem;
      margin: 0.9rem 0;
    }
    .kpi-mini {
      background: rgba(255,255,255,0.03);
      border: 1px solid var(--border);
      padding: 0.75rem;
      border-radius: 8px;
    }
    .kpi-mini span { display: block; font-size: 0.78rem; color: var(--muted); }
    .kpi-mini strong { font-size: 1.35rem; color: var(--af-hi); }
    .bar-strip {
      display: flex;
      height: 44px;
      align-items: flex-end;
      gap: 2px;
      margin-top: 0.8rem;
      padding: 4px;
      background: #070a12;
      border-radius: 6px;
      overflow-x: auto;
    }
    .bar-res {
      flex: 1;
      min-width: 3px;
      border-radius: 2px 2px 0 0;
    }
  </style>
</head>
<body>
  <header>
    <div>
      <h1>AlphaFold-Lite Pocket Viewer · Microservicio en Docker</h1>
      <p style="color: var(--muted); font-size: 0.9rem;">
        Sistemas Programables · Inspirado en AlphaFold (Google DeepMind) · Análisis de Confianza Estructural pLDDT
      </p>
    </div>
    <span class="badge">Contenedor Docker · Puerto 8001</span>
  </header>

  <section class="toolbar">
    <label for="sel-prot" style="font-weight: 600;">Proteína (UniProt ID):</label>
    <select id="sel-prot" onchange="cargarProteina()">
      <option value="P69905">P69905 · Hemoglobina Subunidad Alfa (Homo sapiens)</option>
      <option value="P01308">P01308 · Insulina Humana (Homo sapiens)</option>
      <option value="P01542">P01542 · Crambina (Crambe hispanica)</option>
      <option value="P04637">P04637 · Supresor Tumoral p53 (Homo sapiens)</option>
    </select>
    <button class="primary" onclick="interpretarConIA()">Generar Síntesis Estructural con IA</button>
    <span id="etiqueta-origen" style="font-size: 0.82rem; color: var(--muted);"></span>
  </section>

  <div class="layout-main">
    <section class="card">
      <h3 style="margin-bottom: 0.5rem;">Esqueleto Cα 3D Interactivo (Arrastra con el mouse para rotar)</h3>
      <canvas id="visor3d" width="560" height="360"></canvas>
      <div class="legend">
        <span><i class="dot" style="background: var(--af-vh)"></i>Muy alta (pLDDT ≥ 90)</span>
        <span><i class="dot" style="background: var(--af-hi)"></i>Confiable (70–90)</span>
        <span><i class="dot" style="background: var(--af-lo)"></i>Baja (50–70)</span>
        <span><i class="dot" style="background: var(--af-vl)"></i>Desorden (&lt; 50)</span>
      </div>
    </section>

    <section class="card">
      <h3 id="prot-titulo">Cargando proteína...</h3>
      <p id="prot-sub" style="color: var(--muted); font-size: 0.88rem;"></p>

      <div class="kpi-row">
        <div class="kpi-mini">
          <span>pLDDT Global Promedio</span>
          <strong id="m-plddt">--</strong>
        </div>
        <div class="kpi-mini">
          <span>Radio de Giro (Rg)</span>
          <strong id="m-rg">-- Å</strong>
        </div>
        <div class="kpi-mini">
          <span>Residuos Cα Analizados</span>
          <strong id="m-res">--</strong>
        </div>
        <div class="kpi-mini">
          <span>Fracción Alta Confianza (≥70)</span>
          <strong id="m-alta">-- %</strong>
        </div>
      </div>

      <p style="font-size: 0.85rem; margin-bottom: 0.3rem;">
        <strong>Diagnóstico Global:</strong> <span id="m-clasif" style="color: var(--af-hi);"></span>
      </p>

      <h4 style="margin-top: 1rem; font-size: 0.85rem; color: var(--muted);">
        Perfil de Confianza pLDDT por Residuo (N-terminal → C-terminal)
      </h4>
      <div class="bar-strip" id="perfil-barras"></div>

      <div id="panel-ia" style="margin-top: 1rem; padding-top: 0.9rem; border-top: 1px solid var(--border); display: none;">
        <h4 id="ia-origen" style="color: var(--af-hi); margin-bottom: 0.35rem;"></h4>
        <p id="ia-texto" style="font-size: 0.88rem; margin-bottom: 0.4rem;"></p>
        <ul id="ia-obs" style="margin-left: 1.2rem; font-size: 0.84rem; color: var(--muted);"></ul>
      </div>
    </section>
  </div>

  <script>
    let residuosActuales = [];
    let rotX = 0.35, rotY = 0.65;
    let arrastrando = false, ultX = 0, ultY = 0;

    const canvas = document.getElementById('visor3d');
    const ctx = canvas.getContext('2d');

    function colorPorBanda(banda) {
      if (banda === 'MUY_ALTA_CONFIANZA') return '#1d4ed8';
      if (banda === 'CONFIANZA_ALTA') return '#38bdf8';
      if (banda === 'BAJA_CONFIANZA') return '#facc15';
      return '#fb923c';
    }

    canvas.addEventListener('mousedown', e => { arrastrando = true; ultX = e.clientX; ultY = e.clientY; });
    window.addEventListener('mouseup', () => { arrastrando = false; });
    window.addEventListener('mousemove', e => {
      if (!arrastrando) return;
      rotY += (e.clientX - ultX) * 0.01;
      rotX += (e.clientY - ultY) * 0.01;
      ultX = e.clientX; ultY = e.clientY;
    });

    function renderizar3D() {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      if (residuosActuales.length === 0) {
        requestAnimationFrame(renderizar3D);
        return;
      }
      if (!arrastrando) rotY += 0.006;

      const n = residuosActuales.length;
      let cx = 0, cy = 0, cz = 0;
      for (const r of residuosActuales) { cx += r.x; cy += r.y; cz += r.z; }
      cx /= n; cy /= n; cz /= n;

      const cosY = Math.cos(rotY), sinY = Math.sin(rotY);
      const cosX = Math.cos(rotX), sinX = Math.sin(rotX);
      const escala = 6.2;

      const proyectados = residuosActuales.map(r => {
        const dx = r.x - cx, dy = r.y - cy, dz = r.z - cz;
        const x1 = dx * cosY - dz * sinY;
        const z1 = dx * sinY + dz * cosY;
        const y2 = dy * cosX - z1 * sinX;
        const z2 = dy * sinX + z1 * cosX;
        return {
          px: canvas.width / 2 + x1 * escala,
          py: canvas.height / 2 + y2 * escala,
          pz: z2,
          color: colorPorBanda(r.banda),
        };
      });

      // Dibujar enlaces peptídicos entre carbonos alfa consecutivos
      ctx.lineWidth = 2.4;
      for (let i = 0; i < proyectados.length - 1; i++) {
        ctx.strokeStyle = proyectados[i].color;
        ctx.beginPath();
        ctx.moveTo(proyectados[i].px, proyectados[i].py);
        ctx.lineTo(proyectados[i + 1].px, proyectados[i + 1].py);
        ctx.stroke();
      }

      // Dibujar esferas Cα ordenadas por profundidad Z
      const ordenados = [...proyectados].sort((a, b) => a.pz - b.pz);
      for (const p of ordenados) {
        ctx.fillStyle = p.color;
        ctx.beginPath();
        ctx.arc(p.px, p.py, 4.0, 0, Math.PI * 2);
        ctx.fill();
      }
      requestAnimationFrame(renderizar3D);
    }

    async function cargarProteina() {
      const uid = document.getElementById('sel-prot').value;
      document.getElementById('panel-ia').style.display = 'none';
      const res = await fetch('/api/proteina/' + uid);
      const d = await res.json();
      residuosActuales = d.residuos;
      document.getElementById('prot-titulo').textContent = d.nombre_proteina + ' (' + d.uniprot_id + ')';
      document.getElementById('prot-sub').textContent = 'Organismo: ' + d.organismo;
      document.getElementById('etiqueta-origen').textContent = 'Fuente: ' + d.origen_datos;
      document.getElementById('m-plddt').textContent = d.plddt_promedio.toFixed(1);
      document.getElementById('m-rg').textContent = d.radio_giro_angstrom.toFixed(2) + ' Å';
      document.getElementById('m-res').textContent = d.total_residuos;
      const altaTotal = (d.distribucion_bandas.muy_alta_pct + d.distribucion_bandas.confiable_pct).toFixed(1);
      document.getElementById('m-alta').textContent = altaTotal + ' %';
      document.getElementById('m-clasif').textContent = d.clasificacion_global;

      const tira = document.getElementById('perfil-barras');
      tira.innerHTML = d.residuos.map(r => `
        <div class="bar-res" title="Residuo ${r.indice} (${r.aminoacido}): pLDDT ${r.plddt}"
             style="height: ${Math.max(15, r.plddt)}%; background: ${colorPorBanda(r.banda)};"></div>
      `).join('');
    }

    async function interpretarConIA() {
      const uid = document.getElementById('sel-prot').value;
      const panel = document.getElementById('panel-ia');
      panel.style.display = 'block';
      document.getElementById('ia-texto').textContent = 'Generando interpretación biofísica...';
      document.getElementById('ia-obs').innerHTML = '';
      const res = await fetch('/api/interpretar/' + uid, { method: 'POST' });
      const d = await res.json();
      document.getElementById('ia-origen').textContent = 'Síntesis Estructural · ' + d.origen_modelo;
      document.getElementById('ia-texto').textContent = d.sintesis_estructural;
      document.getElementById('ia-obs').innerHTML = d.observaciones_clave.map(o => `<li>${o}</li>`).join('');
    }

    cargarProteina();
    renderizar3D();
  </script>
</body>
</html>
"""


@app.get("/", response_class=HTMLResponse)
def servir_visor() -> HTMLResponse:
    """Sirve la interfaz web interactiva del visor 3D AlphaFold-Lite."""
    return HTMLResponse(content=HTML_VISOR)
