#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
warmwasser_gui.py
==================
Interaktive GUI zur Berechnung des Wärmeverlusts eines Warmwasserkessels.
Alle Parameter sind per Slider und Eingabefeld konfigurierbar.
Diagramme aktualisieren sich in Echtzeit.

Benötigt:  pip install matplotlib numpy
"""

import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Farben
# ---------------------------------------------------------------------------
C_GREEN  = "#2E7D32"
C_ORANGE = "#E65100"
C_BLUE   = "#1565C0"
C_RED    = "#B71C1C"
C_GREY   = "#90A4AE"
C_BG     = "#F1F8E9"
C_PANEL  = "#E8F5E9"
C_HEADER = "#1B5E20"


# ---------------------------------------------------------------------------
# Physik-Modell
# ---------------------------------------------------------------------------
def compute(V_liter, T_start, T_end, T_amb, n_hours, c_water, rho):
    """Gibt dict mit allen Ergebnissen zurück. None bei ungültigen Eingaben."""
    if T_start <= T_amb or T_end <= T_amb or T_end >= T_start or n_hours <= 0:
        return None
    m        = V_liter * rho
    C_therm  = m * c_water
    n_sec    = n_hours * 3600
    k        = np.log((T_start - T_amb) / (T_end - T_amb)) / n_sec
    UA       = C_therm * k
    Q_J      = C_therm * (T_start - T_end)
    P_avg    = Q_J / n_sec
    P_start  = UA * (T_start - T_amb)
    P_end    = UA * (T_end   - T_amb)
    t_half_h = np.log(2) / k / 3600

    t_s  = np.linspace(0, n_sec, 1000)
    t_h  = t_s / 3600
    T_t  = T_amb + (T_start - T_amb) * np.exp(-k * t_s)
    P_t  = UA * (T_t - T_amb)
    Q_kWh_t = C_therm * (T_start - T_t) / 3_600_000

    return dict(
        m=m, C_therm=C_therm, k=k, k_per_h=k*3600,
        UA=UA, n_sec=n_sec,
        Q_J=Q_J, Q_kJ=Q_J/1000, Q_kWh=Q_J/3_600_000, Q_MJ=Q_J/1e6,
        P_avg=P_avg, P_start=P_start, P_end=P_end,
        t_half_h=t_half_h,
        t_h=t_h, T_t=T_t, P_t=P_t, Q_kWh_t=Q_kWh_t,
        T_start=T_start, T_end=T_end, T_amb=T_amb, n_hours=n_hours,
    )


# ---------------------------------------------------------------------------
# Haupt-App
# ---------------------------------------------------------------------------
class WarmwasserApp:
    PARAMS = [
        # (key, label, unit, from_, to_, resolution, default)
        ("V",       "Volumen",           "L",    10,   2000,  10,   240),
        ("T_start", "Starttemperatur",   "°C",   21,   95,    0.5,  40),
        ("T_end",   "Endtemperatur",     "°C",   18,   90,    0.5,  20),
        ("T_amb",   "Umgebungstemperatur","°C",  -10,  35,    0.5,  17),
        ("n_hours", "Abkühlzeit",        "h",    0.5,  72,    0.5,  8),
        ("c_water", "Spez. Wärmekapazität","J/(kg·K)", 3000, 5000, 10, 4186),
        ("rho",     "Dichte Wasser",     "kg/L",  0.8,  1.1,  0.01, 1.0),
    ]

    def __init__(self, root):
        self.root = root
        self.root.title("Wärmeverlust Warmwasserkessel")
        self.root.configure(bg=C_BG)
        self.root.minsize(1100, 650)

        self._vars = {}      # StringVar für Entry-Felder
        self._sliders = {}   # Scale-Widgets
        self._updating = False

        self._build_ui()
        self._update()

    # ------------------------------------------------------------------ UI --

    def _build_ui(self):
        # ── Top title bar
        title_bar = tk.Frame(self.root, bg=C_HEADER, pady=6)
        title_bar.pack(fill="x")
        tk.Label(title_bar, text="Wärmeverlust Warmwasserkessel",
                 bg=C_HEADER, fg="white",
                 font=("Segoe UI", 14, "bold")).pack(side="left", padx=14)
        tk.Label(title_bar,
                 text="Newtonsches Abkühlungsgesetz  ·  T(t) = T_amb + (T₀ − T_amb) · e^(−k·t)",
                 bg=C_HEADER, fg="#C8E6C9",
                 font=("Segoe UI", 9, "italic")).pack(side="left", padx=8)

        # ── Main horizontal split
        main = tk.Frame(self.root, bg=C_BG)
        main.pack(fill="both", expand=True, padx=6, pady=6)

        # Left control panel
        ctrl = tk.Frame(main, bg=C_PANEL, bd=1, relief="groove", width=310)
        ctrl.pack(side="left", fill="y", padx=(0, 6))
        ctrl.pack_propagate(False)
        self._build_controls(ctrl)

        # Right plot area
        plot_frame = tk.Frame(main, bg=C_BG)
        plot_frame.pack(side="left", fill="both", expand=True)
        self._build_plots(plot_frame)

        # ── Status bar
        self.status_var = tk.StringVar(value="Bereit.")
        tk.Label(self.root, textvariable=self.status_var,
                 bg="#CFD8DC", fg="#37474F",
                 font=("Segoe UI", 9), anchor="w", padx=8
                 ).pack(fill="x", side="bottom")

    def _build_controls(self, parent):
        # ── Section: Parameter
        tk.Label(parent, text="  Parameter", bg=C_HEADER, fg="white",
                 font=("Segoe UI", 10, "bold"), anchor="w"
                 ).pack(fill="x", pady=(0, 4))

        slider_frame = tk.Frame(parent, bg=C_PANEL)
        slider_frame.pack(fill="x", padx=6)

        for key, label, unit, from_, to_, res, default in self.PARAMS:
            self._add_param_row(slider_frame, key, label, unit,
                                from_, to_, res, default)

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8, padx=6)

        # ── Buttons
        btn_frame = tk.Frame(parent, bg=C_PANEL)
        btn_frame.pack(fill="x", padx=6, pady=2)

        tk.Button(btn_frame, text="Zurücksetzen", command=self._reset,
                  bg="#546E7A", fg="white", font=("Segoe UI", 9),
                  relief="flat", padx=8, pady=4
                  ).pack(side="left", padx=(0, 4))
        tk.Button(btn_frame, text="Als PNG speichern", command=self._save_png,
                  bg=C_GREEN, fg="white", font=("Segoe UI", 9, "bold"),
                  relief="flat", padx=8, pady=4
                  ).pack(side="left")

        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=8, padx=6)

        # ── Section: Ergebnisse
        tk.Label(parent, text="  Ergebnisse", bg=C_HEADER, fg="white",
                 font=("Segoe UI", 10, "bold"), anchor="w"
                 ).pack(fill="x")

        res_frame = tk.Frame(parent, bg=C_PANEL)
        res_frame.pack(fill="x", padx=10, pady=6)

        self._result_labels = {}
        results_def = [
            ("Q_kWh",   "Wärmeverlust Q",       "kWh",  C_BLUE,   True),
            ("Q_kJ",    "  = Q",                "kJ",   C_BLUE,   False),
            ("Q_MJ",    "  = Q",                "MJ",   C_BLUE,   False),
            ("P_avg",   "Ø Verlustleistung",    "W",    C_ORANGE, False),
            ("P_start", "Leistung (Anfang)",    "W",    C_ORANGE, False),
            ("P_end",   "Leistung (Ende)",      "W",    C_ORANGE, False),
            ("UA",      "UA-Wert",              "W/K",  C_GREEN,  False),
            ("k_per_h", "Abkühlkonstante k",    "1/h",  C_GREEN,  False),
            ("t_half",  "Halbwertszeit",         "h",    C_GREEN,  False),
        ]
        for key, label, unit, color, bold in results_def:
            row = tk.Frame(res_frame, bg=C_PANEL)
            row.pack(fill="x", pady=1)
            font = ("Segoe UI", 9, "bold") if bold else ("Segoe UI", 9)
            tk.Label(row, text=label, bg=C_PANEL, fg="#333",
                     font=font, width=20, anchor="w").pack(side="left")
            val_lbl = tk.Label(row, text="—", bg=C_PANEL, fg=color,
                               font=("Consolas", 10, "bold"), width=9, anchor="e")
            val_lbl.pack(side="left")
            tk.Label(row, text=unit, bg=C_PANEL, fg="#666",
                     font=("Segoe UI", 9), width=6, anchor="w").pack(side="left")
            self._result_labels[key] = val_lbl

        # Qualitäts-Indikator
        ttk.Separator(parent, orient="horizontal").pack(fill="x", pady=(8,4), padx=6)
        qual_frame = tk.Frame(parent, bg=C_PANEL)
        qual_frame.pack(fill="x", padx=10)
        tk.Label(qual_frame, text="Isoliergüte:", bg=C_PANEL,
                 font=("Segoe UI", 9), fg="#333").pack(side="left")
        self._qual_label = tk.Label(qual_frame, text="—", bg=C_PANEL,
                                    font=("Segoe UI", 9, "bold"), fg=C_GREY)
        self._qual_label.pack(side="left", padx=6)

    def _add_param_row(self, parent, key, label, unit, from_, to_, res, default):
        row = tk.Frame(parent, bg=C_PANEL)
        row.pack(fill="x", pady=3)

        # Label + unit
        tk.Label(row, text=f"{label}", bg=C_PANEL, fg=C_HEADER,
                 font=("Segoe UI", 9, "bold"), width=22, anchor="w"
                 ).grid(row=0, column=0, columnspan=3, sticky="w")

        # Slider
        var = tk.DoubleVar(value=default)
        slider = ttk.Scale(row, from_=from_, to=to_, variable=var,
                           orient="horizontal", length=190,
                           command=lambda v, k=key: self._on_slider(k))
        slider.grid(row=1, column=0, sticky="ew", pady=1)
        self._sliders[key] = slider

        # Entry field
        sv = tk.StringVar(value=str(default))
        entry = ttk.Entry(row, textvariable=sv, width=7,
                          font=("Consolas", 9))
        entry.grid(row=1, column=1, padx=4)
        entry.bind("<Return>",   lambda e, k=key: self._on_entry(k))
        entry.bind("<FocusOut>", lambda e, k=key: self._on_entry(k))

        # Unit label
        tk.Label(row, text=unit, bg=C_PANEL, fg="#666",
                 font=("Segoe UI", 8), width=9, anchor="w"
                 ).grid(row=1, column=2)

        self._vars[key] = {"dvar": var, "svar": sv, "entry": entry,
                           "from_": from_, "to_": to_, "res": res}
        row.columnconfigure(0, weight=1)

    def _build_plots(self, parent):
        self.fig = plt.Figure(figsize=(9, 6.5), facecolor=C_BG)
        self.fig.subplots_adjust(hspace=0.42, wspace=0.32,
                                 left=0.08, right=0.97,
                                 top=0.93, bottom=0.09)
        gs = GridSpec(2, 2, figure=self.fig)
        self.ax_T = self.fig.add_subplot(gs[0, 0])
        self.ax_P = self.fig.add_subplot(gs[0, 1])
        self.ax_Q = self.fig.add_subplot(gs[1, 0])
        self.ax_i = self.fig.add_subplot(gs[1, 1])

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # --------------------------------------------------------- Callbacks ----

    def _on_slider(self, key):
        if self._updating:
            return
        self._updating = True
        v = self._vars[key]["dvar"].get()
        # Round to resolution
        res = self._vars[key]["res"]
        v = round(round(v / res) * res, 6)
        self._vars[key]["svar"].set(f"{v:.4g}")
        self._updating = False
        self._update()

    def _on_entry(self, key):
        if self._updating:
            return
        self._updating = True
        try:
            v = float(self._vars[key]["svar"].get())
            v = max(self._vars[key]["from_"], min(self._vars[key]["to_"], v))
            self._vars[key]["dvar"].set(v)
            self._vars[key]["svar"].set(f"{v:.4g}")
        except ValueError:
            # restore slider value on bad input
            self._vars[key]["svar"].set(
                f"{self._vars[key]['dvar'].get():.4g}")
        self._updating = False
        self._update()

    def _get_params(self):
        return {k: self._vars[k]["dvar"].get() for k in self._vars}

    def _reset(self):
        self._updating = True
        for key, label, unit, from_, to_, res, default in self.PARAMS:
            self._vars[key]["dvar"].set(default)
            self._vars[key]["svar"].set(str(default))
        self._updating = False
        self._update()

    def _save_png(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("Alle Dateien", "*.*")],
            initialfile="warmwasser_verlust.png",
            title="Diagramm speichern")
        if path:
            self.fig.savefig(path, dpi=150, bbox_inches="tight",
                             facecolor=self.fig.get_facecolor())
            self.status_var.set(f"Gespeichert: {path}")

    # --------------------------------------------------------- Compute ------

    def _update(self):
        p = self._get_params()
        r = compute(p["V"], p["T_start"], p["T_end"], p["T_amb"],
                    p["n_hours"], p["c_water"], p["rho"])
        if r is None:
            self.status_var.set(
                "⚠  Ungültige Parameter: T_end und T_start müssen "
                "größer als T_amb sein, und T_start > T_end.")
            for lbl in self._result_labels.values():
                lbl.config(text="—")
            self._qual_label.config(text="—", fg=C_GREY)
            return

        self._update_results(r)
        self._update_plots(r, p)
        self.status_var.set(
            f"  Q = {r['Q_kWh']:.3f} kWh  ·  "
            f"UA = {r['UA']:.2f} W/K  ·  "
            f"Halbwertszeit = {r['t_half_h']:.2f} h  ·  "
            f"Ø Leistung = {r['P_avg']:.1f} W")

    def _update_results(self, r):
        fmts = {
            "Q_kWh":   f"{r['Q_kWh']:.3f}",
            "Q_kJ":    f"{r['Q_kJ']:.1f}",
            "Q_MJ":    f"{r['Q_MJ']:.3f}",
            "P_avg":   f"{r['P_avg']:.1f}",
            "P_start": f"{r['P_start']:.1f}",
            "P_end":   f"{r['P_end']:.1f}",
            "UA":      f"{r['UA']:.2f}",
            "k_per_h": f"{r['k_per_h']:.4f}",
            "t_half":  f"{r['t_half_h']:.2f}",
        }
        for k, v in fmts.items():
            self._result_labels[k].config(text=v)

        # Qualitätsbewertung basierend auf UA-Wert
        ua = r["UA"]
        if ua < 2:
            qual, color = "Sehr gut isoliert", C_GREEN
        elif ua < 5:
            qual, color = "Gut isoliert", "#558B2F"
        elif ua < 15:
            qual, color = "Mässig isoliert", C_ORANGE
        elif ua < 40:
            qual, color = "Schlecht isoliert", "#BF360C"
        else:
            qual, color = "Sehr schlecht isoliert", C_RED
        self._qual_label.config(text=qual, fg=color)

    def _update_plots(self, r, p):
        for ax in (self.ax_T, self.ax_P, self.ax_Q, self.ax_i):
            ax.cla()

        t_h = r["t_h"]
        n_h = r["n_hours"]

        # ---- T(t) ----
        ax = self.ax_T
        ax.plot(t_h, r["T_t"], color=C_GREEN, lw=2.2)
        ax.axhline(r["T_amb"],   color=C_GREY,   lw=1.2, ls="--",
                   label=f"Umgebung {r['T_amb']:.1f} °C")
        ax.axhline(r["T_start"], color=C_ORANGE, lw=1.0, ls=":", alpha=0.8)
        ax.axhline(r["T_end"],   color=C_BLUE,   lw=1.0, ls=":", alpha=0.8)
        ax.fill_between(t_h, r["T_amb"], r["T_t"],
                        alpha=0.15, color=C_GREEN, label="Wärme über Umgebung")
        ax.axvline(r["t_half_h"], color="#9C27B0", lw=1.2, ls="-.",
                   label=f"t½ = {r['t_half_h']:.1f} h")
        ax.annotate(f"{r['T_start']:.1f} °C",
                    xy=(0, r["T_start"]), xytext=(n_h*0.04, r["T_start"]+0.5),
                    color=C_ORANGE, fontsize=8, fontweight="bold")
        ax.annotate(f"{r['T_end']:.1f} °C",
                    xy=(n_h, r["T_end"]),  xytext=(n_h*0.55, r["T_end"]+0.5),
                    color=C_BLUE, fontsize=8, fontweight="bold")
        ax.set_xlabel("Zeit  [h]", fontsize=8)
        ax.set_ylabel("Temperatur  [°C]", fontsize=8)
        ax.set_title("Temperaturverlauf  T(t)", fontweight="bold", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_xlim(0, n_h)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        # ---- P(t) ----
        ax = self.ax_P
        ax.plot(t_h, r["P_t"]/1000, color=C_ORANGE, lw=2.2)
        ax.fill_between(t_h, 0, r["P_t"]/1000, alpha=0.18, color=C_ORANGE)
        ax.axhline(r["P_avg"]/1000, color=C_BLUE, lw=1.4, ls="--",
                   label=f"Ø {r['P_avg']:.0f} W")
        ax.annotate(f"{r['P_start']:.0f} W",
                    xy=(0, r["P_start"]/1000),
                    xytext=(n_h*0.05, r["P_start"]/1000*0.88),
                    color=C_ORANGE, fontsize=8, fontweight="bold")
        ax.annotate(f"{r['P_end']:.0f} W",
                    xy=(n_h, r["P_end"]/1000),
                    xytext=(n_h*0.6, r["P_end"]/1000 + r["P_start"]/1000*0.04),
                    color=C_ORANGE, fontsize=8, fontweight="bold")
        ax.set_xlabel("Zeit  [h]", fontsize=8)
        ax.set_ylabel("Verlustleistung  [kW]", fontsize=8)
        ax.set_title("Momentane Verlustleistung  P(t)", fontweight="bold", fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlim(0, n_h)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        # ---- Q(t) ----
        ax = self.ax_Q
        ax.plot(t_h, r["Q_kWh_t"], color=C_BLUE, lw=2.2)
        ax.fill_between(t_h, 0, r["Q_kWh_t"], alpha=0.15, color=C_BLUE)
        ax.axhline(r["Q_kWh"], color=C_GREEN, lw=1.4, ls="--",
                   label=f"Gesamt {r['Q_kWh']:.3f} kWh")
        ax.annotate(f"{r['Q_kWh']:.3f} kWh\n= {r['Q_kJ']:.0f} kJ",
                    xy=(n_h, r["Q_kWh"]),
                    xytext=(n_h*0.5, r["Q_kWh"]*0.65),
                    color=C_BLUE, fontsize=8, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=C_BLUE, lw=1.0))
        ax.set_xlabel("Zeit  [h]", fontsize=8)
        ax.set_ylabel("Kum. Wärmeverlust  [kWh]", fontsize=8)
        ax.set_title("Kumulierter Wärmeverlust  Q(t)", fontweight="bold", fontsize=9)
        ax.legend(fontsize=7)
        ax.set_xlim(0, n_h)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

        # ---- Info-Panel ----
        ax = self.ax_i
        ax.axis("off")
        ua = r["UA"]
        if ua < 2:     qual, qcol = "Sehr gut isoliert",     C_GREEN
        elif ua < 5:   qual, qcol = "Gut isoliert",          "#558B2F"
        elif ua < 15:  qual, qcol = "Mässig isoliert",       C_ORANGE
        elif ua < 40:  qual, qcol = "Schlecht isoliert",     "#BF360C"
        else:          qual, qcol = "Sehr schlecht isoliert", C_RED

        lines = [
            (C_HEADER, True,  "Eingabe"),
            (C_GREY,   False, f"  Volumen              {p['V']:.0f} L"),
            (C_GREY,   False, f"  T_start / T_end      {p['T_start']:.1f} / {p['T_end']:.1f} °C"),
            (C_GREY,   False, f"  T_Umgebung           {p['T_amb']:.1f} °C"),
            (C_GREY,   False, f"  Abkühlzeit           {p['n_hours']:.1f} h"),
            ("",       False, ""),
            (C_HEADER, True,  "Ergebnisse"),
            (C_BLUE,   True,  f"  Q = {r['Q_kWh']:.3f} kWh  =  {r['Q_kJ']:.0f} kJ"),
            (C_ORANGE, False, f"  Ø Leistung           {r['P_avg']:.1f} W"),
            (C_ORANGE, False, f"  P(0) / P(t_n)        {r['P_start']:.1f} / {r['P_end']:.1f} W"),
            (C_GREEN,  False, f"  UA-Wert              {r['UA']:.2f} W/K"),
            (C_GREEN,  False, f"  k                    {r['k_per_h']:.4f} h⁻¹"),
            (C_GREEN,  False, f"  Halbwertszeit        {r['t_half_h']:.2f} h"),
            ("",       False, ""),
            (qcol,     True,  f"  Isoliergüte: {qual}"),
        ]
        y = 0.97
        for color, bold, text in lines:
            if text == "":
                y -= 0.04
                continue
            font = dict(fontsize=9, fontweight="bold" if bold else "normal")
            ax.text(0.02, y, text, transform=ax.transAxes,
                    color=color, va="top", **font)
            y -= 0.063

        ax.set_title("Zusammenfassung", fontweight="bold", fontsize=9)

        # Highlight box for Q result
        ax.add_patch(plt.Rectangle(
            (0.0, 0.97 - 9*0.063 - 0.01), 1.0, 0.073,
            transform=ax.transAxes,
            facecolor="#E3F2FD", edgecolor="#1565C0", lw=1.2, zorder=0))

        self.fig.suptitle(
            f"Wärmeverlust  ·  {p['V']:.0f} L  ·  "
            f"{p['T_start']:.1f}→{p['T_end']:.1f} °C  ·  "
            f"{p['n_hours']:.1f} h  ·  Umgebung {p['T_amb']:.1f} °C",
            fontsize=10, fontweight="bold", color=C_HEADER, y=0.98)

        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1200x720")
    app = WarmwasserApp(root)
    root.mainloop()
