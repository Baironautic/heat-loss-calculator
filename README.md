# Heat Loss Calculator — Warmwasserkessel

Python tool to calculate and visualise heat loss from a hot water boiler using **Newton's Law of Cooling**.

## Physics Model

```
T(t) = T_amb + (T_start − T_amb) · e^(−k·t)

k  = ln((T_start − T_amb) / (T_end − T_amb)) / t_n     [1/s]
UA = m · c · k                                           [W/K]  thermal conductance
Q  = m · c · (T_start − T_end)                          [J]    total heat loss
t½ = ln(2) / k                                          [s]    half-life
```

## Example

240 L boiler cools from **40 °C → 20 °C** in 8 hours at 17 °C ambient:

| Result | Value |
|---|---|
| Heat loss Q | 20 093 kJ = **5.58 kWh** |
| Avg. power loss | 698 W |
| UA-value | 71 W/K |
| Half-life t½ | 2.72 h |

## Installation (conda)

```bash
conda env create -f environment.yml
conda activate heat-loss-calculator
```

## Usage

### Interactive GUI
```bash
python warmwasser_gui.py
```

All 7 parameters are adjustable via sliders and entry fields. Charts update in real time.

| Parameter | Default | Range |
|---|---|---|
| Volume | 240 L | 10 – 2000 L |
| T_start | 40 °C | 21 – 95 °C |
| T_end | 20 °C | 18 – 90 °C |
| T_ambient | 17 °C | −10 – 35 °C |
| Cooling time | 8 h | 0.5 – 72 h |
| Specific heat capacity | 4186 J/(kg·K) | 3000 – 5000 |
| Water density | 1.0 kg/L | 0.8 – 1.1 |

### Command-line script
```bash
python warmwasser_verlust.py
```
Edit the constants at the top of the file, then run. Saves chart as `warmwasser_verlust.png`.

## Insulation quality (UA-value guide)

| UA [W/K] | Rating |
|---|---|
| < 2 | Very well insulated |
| 2 – 5 | Well insulated |
| 5 – 15 | Moderate |
| 15 – 40 | Poor |
| > 40 | Very poor |

## Files

| File | Description |
|---|---|
| `warmwasser_gui.py` | Interactive tkinter + matplotlib GUI |
| `warmwasser_verlust.py` | Standalone script with static plots |
| `environment.yml` | Conda environment definition |
