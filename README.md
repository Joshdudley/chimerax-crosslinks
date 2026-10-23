# ChimeraX Crosslinks

Display XL-MS crosslinks as CA-to-CA pseudobonds in UCSF ChimeraX. Load intermolecular or intramolecular CSV files, label distances, style dashes, and interactively highlight links that exceed a distance cutoff.

## Setup

Use ChimeraX 1.8 or newer. No additional packages are needed inside ChimeraX. The optional dash-spacing renderer uses internal ChimeraX APIs, so compatibility can vary with version. This is a standalone script, not a ChimeraX Toolshed bundle.

Download this repository and run in the ChimeraX command line:

```text
runscript "/path/to/chimerax-crosslinks/add_crosslinks.py"
```

Load your structures before adding crosslinks. Run the script once per session, and again after reopening a saved session to restore live coloring and optional dash spacing.

## Intermolecular crosslinks

```text
add_crosslinks_inter protein_1='protein_1, model #1/A' protein_2='protein_2, model #2/B' group_name='Inter_protein_1_to_protein_2' data_file='/path/to/intermolecular.csv' crosslink_color='blue' show_distances=true label_decimals=1 xl_cutoff=27
```

The first two CSV columns specify the residue positions in protein_1 and protein_2, respectively. Their headers can be protein names:

```csv
protein_1,protein_2
774,199
780,205
```

## Intramolecular crosslinks

```text
add_crosslinks_intra protein='protein_1, model #1/A' group_name='Intra_protein_1' data_file='/path/to/intramolecular.csv' crosslink_color='blue' show_distances=true label_decimals=1 xl_cutoff=27
```

Intramolecular CSVs require these headers:

```csv
Residue1,Residue2
353,434
450,500
```

`add_crosslinks_intramolecular` is an alias for `add_crosslinks_intra`. Positions may be integers or labels such as `K774`. Extra columns are ignored. The example CSVs in this repository contain synthetic positions to illustrate the format; adapt them to your structures.

## Options shared by both commands

| Option | Default | Meaning |
| --- | --- | --- |
| `group_name` | Required | Crosslink group name. |
| `data_file` | Required | CSV path; quote paths containing spaces. |
| `crosslink_color` | `blue` | Color for links at or below the cutoff. |
| `show_distances` | `false` | Show the nested distance labels. |
| `xl_cutoff` | `27` | Positive CA-to-CA distance cutoff, in angstroms. |
| `color_by_distance` | `true` | Update link colors as structures move. |
| `violation_color` | `red` | Color for distances strictly greater than the cutoff. |
| `radius` | `0.075` | Pseudobond cylinder radius, in angstroms; controls width. |
| `dashes` | `6` | Nonnegative dash count; zero produces solid links. |
| `dash_spacing` | Unset | Approximate dash-plus-gap repeat distance in angstroms; overrides `dashes`. |
| `label_height` | ChimeraX default | Label height in angstroms. |
| `label_size` | ChimeraX default | Positive integer controlling label rendering quality. |
| `label_decimals` | Current session preference | Number of decimal places, including zero; applies to all monitored distances in the session. |
| `output_model` | Automatic | Exact destination model ID, e.g. `'#2.1'`; parent must exist. |
| `replace` | `true` | Replace the existing crosslink group when recreating it. |

With `output_model='#2.1'`, the crosslink group is nested under model #2 and its distance-label child becomes #2.1.1. With replacement enabled, a global pseudobond group at the requested ID is replaced regardless of its name. Protein structures and other model types are protected. Without `output_model`, replacement uses the group name, ignoring capitalization. Global group names must be unique.

Distances exactly equal to `xl_cutoff` retain `crosslink_color`; longer links turn `violation_color` and return to their original color when shortened. This is a visual threshold on straight-line CA distances, not an assessment of crosslink accessibility or conformational populations.

## Calling from a ChimeraX Python script

Python joins adjacent string literals into one command. Keep shared style settings in one place:

```python
from chimerax.core.commands import run as _cmd

label_height = 1.5
label_size = 64
crosslinks_file = "/path/to/intramolecular.csv"

_cmd(session, (
    "add_crosslinks_intra "
    "protein='protein_1, model #1/A' "
    "group_name='Intra_protein_1' "
    f"data_file='{crosslinks_file}' "
    "crosslink_color='blue' show_distances=true "
    f"label_height={label_height} label_size={label_size} "
    "label_decimals=1 dash_spacing=1.0 radius=0.2 xl_cutoff=27"
))
```

This Python example belongs in a Python script run by ChimeraX; the single-line commands above can be pasted directly into its command line.

## Input and display behavior

- CSV numbering must match the loaded structure; no automatic offsets or sequence renumbering are applied.
- Missing or ambiguous residues/CA atoms are skipped and reported. Repeated atom pairs are deduplicated.
- Labels live in one child named `distances`. Hide or show that child in the Model Panel.
- Smaller `dash_spacing` makes dashes and gaps smaller. Counts are rounded to even integers and capped at 2000, so spacing is approximate.
- A `self_dimer` column is ignored. To connect two copies of one protein, use the intermolecular command with two model/chain descriptors and an appropriate CSV.

## Development checks

Outside ChimeraX, install the test dependencies and run:

```bash
python -m pip install -r requirements-dev.txt
python -m unittest discover -s . -p 'test_*.py' -v
```

These tests use a lightweight ChimeraX API double. They check parsing, replacement, chain routing, styling and live color logic; they do not replace testing the actual ChimeraX graphics interface. The script itself uses the NumPy supplied by ChimeraX for array operations.
