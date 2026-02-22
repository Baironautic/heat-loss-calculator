#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
warmwasser_verlust.py
======================
Berechnet und visualisiert den Wärmeverlust eines Warmwasserkessels.

Kessel:               240 L
Abkühlvorgang:        40 °C → 20 °C in n Stunden
Umgebungstemperatur:  17 °C

Physikalisches Modell: Newtonsches Abkühlungsgesetz
  dT/dt = -k · (T - T_amb)
  → T(t) = T_amb + (T_0 - T_amb) · e^(−k·t)

Aus den Randbedingungen T(0)=40, T(n)=20, T_amb=17 folgt:
  k = ln((T_0 − T_amb) / (T_n − T_amb)) / n

Wärmemenge gesamt:  Q = m · c · (T_0 − T_n)
Momentanleistung:   P(t) = m · c · k · (T(t) − T_amb)
Thermischer Leitwert: UA = P(t) / (T(t) − T_amb) = m · c · k  [W/K]
"""

import sys
import numpy as np
# Force UTF-8 on Windows console
if sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch

# ---------------------------------------------------------------------------
# Parameter (hier anpassen)
# ---------------------------------------------------------------------------
V_liter      = 240        # Kesselvolumen [L]
T_start      = 40.0       # Starttemperatur [°C]
T_end        = 20.0       # Endtemperatur nach n Stunden [°C]
T_amb        = 17.0       # Umgebungstemperatur [°C]
n_hours      = 8.0        # Abkühlzeit [h]   ← hier ändern

# ---------------------------------------------------------------------------
# Physikalische Konstanten
# ---------------------------------------------------------------------------
rho_water    = 1.0        # Dichte Wasser [kg/L]
c_water      = 4186.0     # Spez. Wärmekapazität Wasser [J/(kg·K)]

# ---------------------------------------------------------------------------
# Berechnungen
# ---------------------------------------------------------------------------
m            = V_liter * rho_water                   # Masse [kg]
C_therm      = m * c_water                           # Thermische Kapazität [J/K]

# Abkühlkonstante k aus T(n) = T_amb + (T_0-T_amb)*exp(-k*n)
n_sec        = n_hours * 3600                        # [s]
k            = np.log((T_start - T_amb) /
                      (T_end   - T_amb)) / n_sec     # [1/s]

# Thermischer Leitwert (UA-Wert)
UA           = C_therm * k                           # [W/K]

# Gesamte abgegebene Wärmemenge
Q_total_J    = C_therm * (T_start - T_end)           # [J]
Q_total_kWh  = Q_total_J / 3_600_000                # [kWh]
Q_total_kJ   = Q_total_J / 1000                     # [kJ]

# Durchschnittliche Verlustleistung
P_avg        = Q_total_J / n_sec                     # [W]

# Momentanleistung bei Start und Ende
P_start      = UA * (T_start - T_amb)               # [W]
P_end        = UA * (T_end   - T_amb)               # [W]

# Zeitachse
t_s          = np.linspace(0, n_sec, 2000)          # [s]
t_h          = t_s / 3600                            # [h]

# Temperaturverlauf T(t)
T_t          = T_amb + (T_start - T_amb) * np.exp(-k * t_s)

# Momentane Verlustleistung P(t) = UA · (T(t) - T_amb)
P_t          = UA * (T_t - T_amb)                   # [W]

# Kumulierter Wärmeverlust Q(t) = C_therm · (T_0 − T(t))
Q_t_J        = C_therm * (T_start - T_t)            # [J]
Q_t_kWh      = Q_t_J / 3_600_000                   # [kWh]

# Halbwertszeit der Abkühlung (T(t½) = T_amb + (T_0-T_amb)/2)
t_half_s     = np.log(2) / k
t_half_h     = t_half_s / 3600

# ---------------------------------------------------------------------------
# Konsolenausgabe
# ---------------------------------------------------------------------------
print("=" * 55)
print("  WÄRMEVERLUST WARMWASSERKESSEL")
print("=" * 55)
print(f"  Volumen           : {V_liter} L")
print(f"  Wassermasse       : {m:.0f} kg")
print(f"  Starttemperatur   : {T_start:.1f} °C")
print(f"  Endtemperatur     : {T_end:.1f} °C  (nach {n_hours:.1f} h)")
print(f"  Umgebung          : {T_amb:.1f} °C")
print("-" * 55)
print(f"  Abkuehlkonstante k: {k*3600:.4f}  1/h  ({k:.6f} 1/s)")
print(f"  UA-Wert (Kessel)  : {UA:.2f} W/K")
print(f"  Halbwertszeit     : {t_half_h:.2f} h")
print("-" * 55)
print(f"  Gesamtwärmeverlust: {Q_total_kJ:.1f} kJ")
print(f"                      {Q_total_kWh:.3f} kWh")
print(f"                      {Q_total_J/1e6:.3f} MJ")
print(f"  Ø Verlustleistung : {P_avg:.1f} W")
print(f"  Leistung (Anfang) : {P_start:.1f} W")
print(f"  Leistung (Ende)   : {P_end:.1f} W")
print("=" * 55)

# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 9), facecolor="#F8F9FA")
fig.suptitle(
    f"Wärmeverlust Warmwasserkessel  ·  {V_liter} L  ·  {T_start:.0f}→{T_end:.0f} °C "
    f"in {n_hours:.0f} h  ·  Umgebung {T_amb:.0f} °C",
    fontsize=13, fontweight="bold", color="#1B5E20", y=0.98
)

gs = gridspec.GridSpec(2, 2, figure=fig,
                       hspace=0.42, wspace=0.35,
                       left=0.08, right=0.97, top=0.93, bottom=0.08)

ax_T  = fig.add_subplot(gs[0, 0])   # Temperaturverlauf
ax_P  = fig.add_subplot(gs[0, 1])   # Verlustleistung
ax_Q  = fig.add_subplot(gs[1, 0])   # Kumulierter Wärmeverlust
ax_i  = fig.add_subplot(gs[1, 1])   # Info-Box

GREEN  = "#2E7D32"
ORANGE = "#E65100"
BLUE   = "#1565C0"
GREY   = "#90A4AE"

# ---- Panel 1: Temperaturverlauf ----
ax_T.plot(t_h, T_t, color=GREEN, lw=2.5, label="T(t)")
ax_T.axhline(T_amb, color=GREY,   lw=1.2, ls="--", label=f"Umgebung {T_amb:.0f} °C")
ax_T.axhline(T_start, color=ORANGE, lw=1.0, ls=":", alpha=0.7)
ax_T.axhline(T_end,   color=BLUE,   lw=1.0, ls=":", alpha=0.7)
ax_T.fill_between(t_h, T_amb, T_t, alpha=0.15, color=GREEN, label="Nutzbare Wärme über Umgebung")
ax_T.axvline(t_half_h, color="#9C27B0", lw=1.2, ls="-.",
             label=f"Halbwertszeit {t_half_h:.1f} h")

ax_T.annotate(f"{T_start:.0f} °C",  xy=(0, T_start), xytext=(0.15, T_start+0.6),
              color=ORANGE, fontsize=9, fontweight="bold")
ax_T.annotate(f"{T_end:.0f} °C",    xy=(n_hours, T_end), xytext=(n_hours*0.6, T_end+0.6),
              color=BLUE, fontsize=9, fontweight="bold")
ax_T.annotate(f"{T_amb:.0f} °C",    xy=(n_hours, T_amb), xytext=(n_hours*0.6, T_amb+0.6),
              color=GREY, fontsize=9)

ax_T.set_xlabel("Zeit  [h]")
ax_T.set_ylabel("Temperatur  [°C]")
ax_T.set_title("Temperaturverlauf  T(t)", fontweight="bold")
ax_T.legend(fontsize=8, loc="upper right")
ax_T.set_xlim(0, n_hours)
ax_T.grid(True, alpha=0.3)

# ---- Panel 2: Verlustleistung ----
ax_P.plot(t_h, P_t, color=ORANGE, lw=2.5)
ax_P.fill_between(t_h, 0, P_t, alpha=0.2, color=ORANGE)
ax_P.axhline(P_avg, color=BLUE, lw=1.5, ls="--",
             label=f"Ø Leistung  {P_avg:.1f} W")

ax_P.annotate(f"P(0) = {P_start:.1f} W", xy=(0, P_start),
              xytext=(n_hours*0.1, P_start*0.88),
              color=ORANGE, fontsize=9, fontweight="bold")
ax_P.annotate(f"P({n_hours:.0f}h) = {P_end:.1f} W", xy=(n_hours, P_end),
              xytext=(n_hours*0.55, P_end + P_start*0.05),
              color=ORANGE, fontsize=9, fontweight="bold")

ax_P.set_xlabel("Zeit  [h]")
ax_P.set_ylabel("Verlustleistung  [W]")
ax_P.set_title("Momentane Verlustleistung  P(t) = UA · (T−T_amb)", fontweight="bold")
ax_P.legend(fontsize=9)
ax_P.set_xlim(0, n_hours)
ax_P.set_ylim(bottom=0)
ax_P.grid(True, alpha=0.3)

# ---- Panel 3: Kumulierter Wärmeverlust ----
ax_Q.plot(t_h, Q_t_kWh, color=BLUE, lw=2.5)
ax_Q.fill_between(t_h, 0, Q_t_kWh, alpha=0.15, color=BLUE)
ax_Q.axhline(Q_total_kWh, color=GREEN, lw=1.5, ls="--",
             label=f"Gesamt {Q_total_kWh:.3f} kWh")

ax_Q.annotate(f"{Q_total_kWh:.3f} kWh\n= {Q_total_kJ:.0f} kJ",
              xy=(n_hours, Q_total_kWh),
              xytext=(n_hours*0.55, Q_total_kWh*0.7),
              color=BLUE, fontsize=9, fontweight="bold",
              arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))

ax_Q.set_xlabel("Zeit  [h]")
ax_Q.set_ylabel("Kumulierter Verlust  [kWh]")
ax_Q.set_title("Kumulierter Wärmeverlust  Q(t)", fontweight="bold")
ax_Q.legend(fontsize=9)
ax_Q.set_xlim(0, n_hours)
ax_Q.set_ylim(bottom=0)
ax_Q.grid(True, alpha=0.3)

# ---- Panel 4: Ergebnis-Zusammenfassung ----
ax_i.axis("off")
summary_lines = [
    ("Gegebene Größen", None, True),
    (f"  Volumen",              f"{V_liter} L",           False),
    (f"  T_start",              f"{T_start:.1f} °C",      False),
    (f"  T_end",                f"{T_end:.1f} °C",        False),
    (f"  T_Umgebung",           f"{T_amb:.1f} °C",        False),
    (f"  Abkühlzeit",           f"{n_hours:.1f} h",       False),
    ("", None, False),
    ("Berechnete Größen", None, True),
    (f"  Abkühlkonstante k",    f"{k*3600:.4f} h⁻¹",     False),
    (f"  UA-Wert",              f"{UA:.2f} W/K",          False),
    (f"  Halbwertszeit",        f"{t_half_h:.2f} h",      False),
    ("", None, False),
    ("Wärmeverlust", None, True),
    (f"  Q gesamt",             f"{Q_total_kJ:.1f} kJ",   False),
    (f"  Q gesamt",             f"{Q_total_kWh:.3f} kWh", False),
    (f"  Q gesamt",             f"{Q_total_J/1e6:.3f} MJ",False),
    (f"  Ø Verlustleistung",    f"{P_avg:.1f} W",         False),
    (f"  P bei Start",          f"{P_start:.1f} W",       False),
    (f"  P bei Ende",           f"{P_end:.1f} W",         False),
]

y = 0.97
for label, value, is_header in summary_lines:
    if is_header:
        ax_i.text(0.02, y, label, transform=ax_i.transAxes,
                  fontsize=10, fontweight="bold", color=GREEN,
                  va="top")
    elif label == "":
        pass
    else:
        ax_i.text(0.03, y, label, transform=ax_i.transAxes,
                  fontsize=9, color="#333333", va="top")
        ax_i.text(0.72, y, value, transform=ax_i.transAxes,
                  fontsize=9, color="#1565C0", va="top", fontweight="bold",
                  ha="right")
    y -= 0.055

# Highlight box around the main result
highlight_y = 0.97 - 12 * 0.055
ax_i.add_patch(plt.Rectangle((0.0, highlight_y - 0.005),
                               1.0, 3 * 0.055 + 0.01,
                               transform=ax_i.transAxes,
                               facecolor="#E8F5E9", edgecolor="#4CAF50",
                               lw=1.5, zorder=0))

ax_i.set_title("Ergebniszusammenfassung", fontweight="bold")

plt.savefig("C:/claude/warmwasser_verlust.png", dpi=150, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print("\nDiagramm gespeichert: C:/claude/warmwasser_verlust.png")
plt.show()
