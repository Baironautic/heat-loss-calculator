# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python heat loss calculator for a hot water boiler (Warmwasserkessel). Models cooling using Newton's law of cooling and visualises results with matplotlib.

## Running the Scripts

```bash
# Interactive GUI (recommended)
python warmwasser_gui.py

# Command-line calculation + static chart
python warmwasser_verlust.py
```

## Dependencies

No `requirements.txt` or `environment.yml` yet. Install manually:

```bash
pip install matplotlib numpy
```

## Architecture

### warmwasser_verlust.py — standalone script

Single-file calculation + plot. Parameters are hardcoded constants at the top of the file (edit and re-run to change values). Saves output to `warmwasser_verlust.png`.

Key constants to adjust:
```python
V_liter   = 240    # Kesselvolumen [L]
T_start   = 40.0   # Starttemperatur [°C]
T_end     = 20.0   # Endtemperatur [°C]
T_amb     = 17.0   # Umgebungstemperatur [°C]
n_hours   = 8.0    # Abkühlzeit [h]
```

### warmwasser_gui.py — interactive GUI

Two-class design:

- `compute(V, T_start, T_end, T_amb, n_hours, c_water, rho)` — pure physics function, no GUI imports. Returns a dict of all results and time-series arrays. Returns `None` for invalid parameter combinations.
- `WarmwasserApp` — tkinter + embedded matplotlib. All 7 parameters are exposed as sliders + Entry fields. Plots re-render on every parameter change via `_update()` → `compute()` → `_update_plots()`.

### Physical model

Newton's law of cooling:
```
T(t) = T_amb + (T_start - T_amb) · exp(-k·t)
k    = ln((T_start - T_amb) / (T_end - T_amb)) / n_sec
UA   = m · c · k          [W/K]  — thermal conductance of the vessel
t½   = ln(2) / k          [s]    — half-life of temperature difference
Q    = m · c · (T_start - T_end) [J]
```

### Known model limitation

`k` (and therefore t½ and UA) is derived from the fixed measured pair (T_end, n_hours). Changing T_amb while keeping T_end fixed alters the computed k — making the boiler appear better or worse insulated depending on ambient temperature. This is a modelling artefact of using T_end as a fixed input rather than UA.
