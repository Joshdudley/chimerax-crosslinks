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

## Interactive coloring for existing pseudobond groups

Use `color_crosslinks` to apply live distance coloring to an existing pseudobond group created by any method. Load `add_crosslinks.py` with `runscript` first, then run:

```text
color_crosslinks model='#3' xl_cutoff=27 violation_color='red'
```

Use the exact pseudobond group ID, such as `#3` or `#2.1`, rather than a protein model, a parent container, or the group's distance-label child. This command does not require a CSV and does not recreate bonds or labels.

By default, each bond retains its original color at or below the cutoff and turns red strictly above it. Colors update as the endpoint atoms or their models move. For groups with atom-colored bond halves, the original half-color behavior is preserved below the cutoff; the whole bond uses the violation color above it. Distances are measured between the actual endpoint atoms, which need not be CA atoms.

To use one color for all bonds below the cutoff:

```text
color_crosslinks model='#3' xl_cutoff=27 crosslink_color='blue' violation_color='red'
```

To disable live coloring and restore the colors and half-color settings saved when standalone coloring was first enabled:

```text
color_crosslinks model='#3' enabled=false
```

| Option | Default | Meaning |
| --- | --- | --- |
| `model` | Required | One exact pseudobond group model ID. |
| `xl_cutoff` | `27` | Positive finite distance cutoff in angstroms. |
| `violation_color` | `red` | Color used strictly above the cutoff. |
| `crosslink_color` | Unset | Optional uniform color at or below the cutoff; omit to use original individual colors. |
| `enabled` | `true` | Enable coloring, or disable and restore saved original styles. |

Repeated calls can change the cutoff or colors without losing the saved originals. Bonds added to the group while coloring is enabled are included on the next graphics update; deleted bonds are no longer tracked. Manual recoloring is overridden while live coloring is enabled. Disable first if you want to establish new original colors. Settings and saved original styles are stored as session attributes; reload this script when reopening a saved session to reactivate updates.

The existing `add_crosslinks_inter` and `add_crosslinks_intra` commands retain their live-coloring options and share the same update mechanism. You can also apply `color_crosslinks` to their groups; its saved originals are the bond styles present when the standalone command is first enabled.

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

These tests use a lightweight ChimeraX API double. They check parsing, replacement, chain routing, styling, live color logic, and standalone color restoration; they do not replace testing the actual ChimeraX graphics interface. The script itself uses the NumPy supplied by ChimeraX for array operations.
