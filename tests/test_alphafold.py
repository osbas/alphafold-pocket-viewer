"""Pruebas unitarias deterministas para AlphaFold-Lite Pocket Viewer."""

from __future__ import annotations

from app.main import (
    analizar_estructura_pdb,
    clasificar_banda_plddt,
    generar_interpretacion_ia,
    generar_pdb_plan_b,
)


def ejecutar_pruebas() -> None:
    """Verifica la clasificación de bandas pLDDT, síntesis PDB y cálculo de Rg."""
    assert clasificar_banda_plddt(95.0) == "MUY_ALTA_CONFIANZA"
    assert clasificar_banda_plddt(82.0) == "CONFIANZA_ALTA"
    assert clasificar_banda_plddt(61.0) == "BAJA_CONFIANZA"
    assert clasificar_banda_plddt(41.0) == "DESORDEN_ESTRUCTURAL"

    pdb_hemoglobina = generar_pdb_plan_b("P69905")
    assert "ATOM" in pdb_hemoglobina
    reporte = analizar_estructura_pdb("P69905", pdb_hemoglobina, "PLAN_B_TEST")
    assert reporte.total_residuos == 142
    assert 80.0 <= reporte.plddt_promedio <= 99.0
    assert reporte.radio_giro_angstrom > 5.0

    suma_pct = (
        reporte.distribucion_bandas.muy_alta_pct
        + reporte.distribucion_bandas.confiable_pct
        + reporte.distribucion_bandas.baja_pct
        + reporte.distribucion_bandas.muy_baja_pct
    )
    assert abs(suma_pct - 100.0) <= 0.5

    interp = generar_interpretacion_ia(reporte)
    assert len(interp.observaciones_clave) == 3
    print("Todas las pruebas unitarias de AlphaFold-Lite finalizaron con éxito.")


if __name__ == "__main__":
    ejecutar_pruebas()
