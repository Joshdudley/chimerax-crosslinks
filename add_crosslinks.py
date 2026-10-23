"""Load this file in ChimeraX with: runscript /path/to/add_crosslinks.py

Commands use quoted key=value options, for example:
add_crosslinks_inter protein_1='protein_1, model #1.1/A' protein_2='protein_2, model #1.2/B' group_name='protein_1_to_protein_2' data_file='/path/intermolecular.csv' crosslink_color='blue' show_distances=false
add_crosslinks_intra protein='protein_1, model #1.1/A' group_name='protein_1_internal' data_file='/path/intramolecular_protein_1.csv' crosslink_color='blue' show_distances=false

CSV positions must already match the loaded structure: no automatic offset or
sequence renumbering is applied. Distances are straight-line CA-to-CA distances.
The CSV's first two columns provide endpoints; protein labels describe the
chosen chains and need not equal CSV headers. Extra columns, including
self_dimer, are ignored. Missing CA atoms/residues are skipped and reported.
If data_file has no extension, a matching path with .csv appended is accepted.

Pseudobond models have a child named 'distances'. It contains updating distance
labels even when show_distances=false hides that child. The links remain visible.
Toggle those children in the Model Panel to show/hide labels later.

To draw links between two copies of the same protein, use add_crosslinks_inter
with two different model/chain descriptors and the appropriate CSV. The inter
command accepts Residue1/Residue2 CSVs as well as protein-named columns.

output_model='#2.1' gives the crosslink group that exact model ID and nests it
under existing model #2. Its distance labels become #2.1.1. With replace=true
(the default), an existing global pseudobond group at #2.1 is deleted together
with its distance labels, regardless of name, and recreated using group_name.
Protein structures and other model types cannot be replaced this way.
The parent must already exist. Input and destination validation happen before
replacement. A group_name already used at another ID must be changed because
global pseudobond names are unique. Omit output_model to replace by group_name,
ignoring capitalization, and use automatic top-level numbering.
Optional options for either command: replace=false (refuse replacement),
radius=0.075, dashes=6. Defaults: blue, no labels, replace=true.
label_height=1.5 (Angstroms) and label_size=64 (rendering quality) optionally style the
nested distance labels; omitted values keep ChimeraX's current label defaults.
label_decimals=1 displays one digit after the decimal point (7.2Å).
Zero is allowed; omission keeps the current preference.
This uses ChimeraX's session-wide distance precision, so existing monitored
distances also adopt it and retain it as structures move. Coordinates and
computed distances are not rounded.
Distance-based coloring is on by default for both commands:
color_by_distance=true xl_cutoff=27 violation_color='red'
Links strictly longer than 27 Angstroms turn red; at or below the cutoff they
return to crosslink_color. This updates interactively during protein/model
movement, works with fixed dashes and dash_spacing, and uses unrounded CA-to-CA
distances. Set color_by_distance=false to keep one color for the whole group.
After reopening a saved session, load this script to reactivate live coloring.
The Log reports the initial minimum, mean and maximum CA-to-CA distances.
dash_spacing=1.0 optionally targets a 1-Angstrom dash-plus-gap repeat distance,
overriding dashes. Smaller spacing creates shorter dashes and gaps. Counts are
rounded to even integers, so spacing is approximate, especially for short links.
Spacing mode uses internal graphics drawings, keeping just one pseudobond model
with one distances child in the Model Panel. Endpoints, labels and dash spacing
update when structures move. After reopening a saved session, run this script
again to restore the spacing renderer (or load it before opening the session).
The renderer uses ChimeraX internal graphics APIs; compatibility may depend on
ChimeraX version.
Re-run this script once per ChimeraX session to register commands.
Requires ChimeraX 1.8+ (the pbond command). No extra Python packages required.
"""

import csv
import math
import re
import shlex
from pathlib import Path


class CrosslinkInputError(ValueError):
    """A malformed CSV row or command option."""


_PROTEIN_SPEC = re.compile(
    r"^\s*(?P<label>[^,]+),\s*model\s+#(?P<model>\d+(?:\.\d+)*)/"
    r"(?P<chain>[A-Za-z0-9_-]+)\s*$", re.IGNORECASE
)
_POSITION = re.compile(r"^(?:[A-Za-z])?(\d+)$")


def parse_protein(value):
    match = _PROTEIN_SPEC.fullmatch(value)
    if match is None:
        raise CrosslinkInputError(
            "Use a protein descriptor such as 'protein_1, model #1.1/A'."
        )
    return match['label'].strip(), match['model'], match['chain']


def parse_bool(value):
    if isinstance(value, bool):
        return value
    token = str(value).strip().lower()
    if token in ('true', 't', '1'):
        return True
    if token in ('false', 'f', '0'):
        return False
    raise CrosslinkInputError("Boolean values must be true or false.")


def parse_position(value, line):
    match = _POSITION.fullmatch(str(value).strip())
    if match is None or int(match[1]) < 1:
        raise CrosslinkInputError(
            f"CSV line {line}: expected a positive residue number or label such as K774; got {value!r}."
        )
    return int(match[1])


def read_pairs(data_file, intramolecular=False):
    """Return endpoint headers and (CSV line, residue1, residue2) rows."""
    path = Path(data_file).expanduser()
    if not path.exists() and not path.suffix:
        path = path.with_suffix('.csv')
    try:
        with path.open(newline='', encoding='utf-8-sig') as stream:
            reader = csv.DictReader(stream)
            headers = reader.fieldnames or []
            if len(headers) < 2 or len(set(headers)) != len(headers):
                raise CrosslinkInputError("CSV needs at least two distinct column headers.")
            if intramolecular:
                columns = ('Residue1', 'Residue2')
                if any(name not in headers for name in columns):
                    raise CrosslinkInputError("Intramolecular CSV needs Residue1 and Residue2 columns.")
            else:
                columns = tuple(headers[:2])
                if 'self_dimer' in columns:
                    raise CrosslinkInputError("The first two CSV columns must contain residue endpoints.")
            pairs = []
            for row in reader:
                line = reader.line_num
                if None in row or any(row.get(column) is None for column in columns):
                    raise CrosslinkInputError(f"Malformed CSV line {line}.")
                pairs.append((line, parse_position(row[columns[0]], line),
                              parse_position(row[columns[1]], line)))
    except (OSError, UnicodeError, csv.Error) as error:
        raise CrosslinkInputError(f"Could not read {path}: {error}") from error
    if not pairs:
        raise CrosslinkInputError(f"No crosslink rows found in {path}.")
    return columns, pairs


def parse_options(text, intramolecular=False):
    try:
        tokens = shlex.split(text, comments=False, posix=True)
    except ValueError as error:
        raise CrosslinkInputError(str(error)) from error
    allowed = {'group_name', 'data_file', 'crosslink_color', 'show_distances',
               'replace', 'radius', 'dashes', 'dash_spacing', 'label_height', 'label_size', 'output_model', 'label_decimals', 'color_by_distance', 'xl_cutoff', 'violation_color'}
    required = {'group_name', 'data_file'}
    if intramolecular:
        allowed.add('protein')
        required.add('protein')
    else:
        allowed.update(('protein_1', 'protein_2'))
        required.update(('protein_1', 'protein_2'))
    values = {}
    for token in tokens:
        key, equals, value = token.partition('=')
        if not equals or key not in allowed:
            raise CrosslinkInputError(f"Unknown option {token!r}. Use key=value with quoted strings.")
        if key in values or not value.strip():
            raise CrosslinkInputError(f"Option {key!r} is repeated or empty.")
        values[key] = value
    missing = required - values.keys()
    if missing:
        raise CrosslinkInputError('Missing options: ' + ', '.join(sorted(missing)))
    for key in ('show_distances', 'replace', 'color_by_distance'):
        if key in values:
            values[key] = parse_bool(values[key])
    try:
        for key in ('radius', 'label_height', 'dash_spacing', 'xl_cutoff'):
            if key in values:
                values[key] = float(values[key])
                if not math.isfinite(values[key]) or values[key] <= 0:
                    raise ValueError
        if 'label_decimals' in values:
            values['label_decimals'] = int(values['label_decimals'])
            if values['label_decimals'] < 0:
                raise ValueError
        if 'dashes' in values:
            values['dashes'] = int(values['dashes'])
            if values['dashes'] < 0:
                raise ValueError
        if 'label_size' in values:
            values['label_size'] = int(values['label_size'])
            if values['label_size'] <= 0:
                raise ValueError
    except ValueError as error:
        raise CrosslinkInputError('radius/label_height/dash_spacing/xl_cutoff must be positive and finite; '
                                 'label_size must be a positive integer; dashes and label_decimals must be nonnegative integers.') from error
    return values


def _chain_residues(session, descriptor):
    from chimerax.atomic import AtomicStructure
    label, model_id, chain_id = parse_protein(descriptor)
    structures = [model for model in session.models.list(type=AtomicStructure)
                  if model.id_string == model_id]
    if len(structures) != 1:
        raise CrosslinkInputError(f"Could not find exactly one atomic model #{model_id} for {label}.")
    residues = [r for r in structures[0].residues if r.chain_id == chain_id]
    if not residues:
        raise CrosslinkInputError(f"Chain #{model_id}/{chain_id} was not found for {label}.")
    index = {}
    for residue in residues:
        index.setdefault(residue.number, []).append(residue)
    return label, f'#{model_id}/{chain_id}', index


def _ca_atom(chain, position):
    residues = chain[2].get(position, [])
    if len(residues) != 1:
        reason = 'missing residue' if not residues else 'ambiguous insertion-coded residue'
        return None, f'{chain[1]}:{position}: {reason}'
    atom = residues[0].find_atom('CA')
    return atom, None if atom is not None else f'{chain[1]}:{position}: missing CA'


def dashes_for_spacing(length, spacing):
    """Approximate length/count = requested repeat spacing with even counts."""
    ratio = length / spacing
    if not math.isfinite(ratio) or ratio > 2000:
        raise CrosslinkInputError('dash_spacing is too small for these links; increase it to use at most 2000 dashes per link.')
    return max(2, 2 * int(ratio / 2 + 0.5))


class _SpacingRenderer:
    """Draw count-specific cylinders without creating any additional Models.

    Native pseudobonds remain intact for distance monitors and atom picking.
    This adapter intentionally isolates the internal ChimeraX graphics API.
    """

    def __init__(self, group, spacing):
        self.group, self.spacing = group, spacing
        self.drawings = {}
        self.original_update = group._update_graphics
        self.original_planes_pick = group._pseudobonds_planes_pick

    def update(self, changes=None):
        group = self.group
        if group.deleted:
            return
        self.original_update(group._ALL_CHANGE if changes is None else changes)
        native = group._pbond_drawing
        if native is None:
            self.clear()
            return
        # Hide the native cylinders only, leaving bonds and their labels shown.
        native.display = False
        bonds = native.visible_bonds
        from numpy import array, linalg
        atoms1, atoms2 = bonds.atoms
        lengths = linalg.norm(atoms1.pb_scene_coords - atoms2.pb_scene_coords, axis=1)
        # Cap after movement, so an extreme interactive displacement never
        # raises from a graphics trigger or allocates unbounded geometry.
        counts = array([dashes_for_spacing(min(float(n) / self.spacing, 2000), 1)
                        for n in lengths])
        active = set(int(n) for n in counts)
        for count in set(self.drawings) - active:
            group.remove_drawing(self.drawings.pop(count))
        from chimerax.atomic.structure import (BondsDrawing, PickedPseudobond,
                                              PickedPseudobonds, _selected_bond_cylinders)
        from chimerax.surface import dashed_cylinder_geometry
        for count in active:
            drawing = self.drawings.get(count)
            if drawing is None:
                drawing = BondsDrawing('crosslink dashes', PickedPseudobond, PickedPseudobonds)
                self.drawings[count] = drawing
                group.add_drawing(drawing)
            sides = group._cylinder_sides
            if getattr(drawing, '_xlms_sides', None) != sides:
                drawing.set_geometry(*dashed_cylinder_geometry(count // 2, height=0.5, nc=sides))
                drawing._xlms_sides = sides
            subset = bonds.filter(counts == count)
            drawing.visible_bonds = subset
            drawing._visible_atoms = subset.atoms
            drawing.positions = group._update_positions(subset, subset.atoms)
            drawing.colors = subset.half_colors
            drawing.highlighted_positions = _selected_bond_cylinders(subset)

    def planes_pick(self, planes):
        # The native implementation picks through _pbond_drawing. Temporarily
        # point it at each visible drawing so rectangle selection still works.
        group = self.group
        native = group._pbond_drawing
        picks = []
        try:
            for drawing in self.drawings.values():
                group._pbond_drawing = drawing
                picks.extend(self.original_planes_pick(planes))
        finally:
            group._pbond_drawing = native
        return picks

    def clear(self):
        for drawing in self.drawings.values():
            self.group.remove_drawing(drawing)
        self.drawings.clear()


def _install_spacing_renderer(group, spacing):
    # Save only the scalar spacing attribute; Drawings and Python hooks are
    # rebuilt when this script is loaded after restoring a saved session.
    group.register_attr(group.session, 'xlms_dash_spacing', 'XLMS crosslinks', attr_type=float)
    group.xlms_dash_spacing = float(spacing)
    if getattr(group, '_xlms_renderer', None) is not None:
        group._xlms_renderer.spacing = float(spacing)
        group._xlms_renderer.update()
        return
    renderer = _SpacingRenderer(group, spacing)
    group._xlms_renderer = renderer
    group._update_graphics = renderer.update
    group._pseudobonds_planes_pick = renderer.planes_pick
    renderer.update()
    group.redraw_needed(shape_changed=True)


def _restore_spacing_renderers(session):
    from chimerax.atomic import PseudobondGroup
    for group in session.models.list(type=PseudobondGroup):
        spacing = getattr(group, 'xlms_dash_spacing', None)
        if spacing is not None:
            _install_spacing_renderer(group, spacing)
    _ensure_distance_color_updates(session)
    _update_distance_colors(session)


def _parse_color_rgba(session, value, option):
    from chimerax.core.commands import ColorArg, StringArg
    color, _, remaining = ColorArg.parse(StringArg.unparse(value), session)
    if remaining.strip():
        raise CrosslinkInputError(f'Specify a single {option}.')
    return tuple(int(channel) for channel in color.uint8x4())


def _update_group_distance_colors(group):
    """Apply the cutoff to true CA distances, including parent model motion."""
    if group.deleted or not getattr(group, 'xlms_color_by_distance', False):
        return
    bonds = group.pseudobonds
    if len(bonds) == 0:
        return
    import numpy as np
    base = np.frombuffer(bytes.fromhex(group.xlms_base_rgba), dtype=np.uint8)
    red = np.frombuffer(bytes.fromhex(group.xlms_violation_rgba), dtype=np.uint8)
    colors = np.tile(base, (len(bonds), 1))
    # Use unrounded distances, not the formatted label or ribbon attachment.
    colors[bonds.lengths > group.xlms_distance_cutoff] = red
    halfbond = bool(np.any(bonds.halfbonds))
    if halfbond:
        bonds.halfbonds = False
    if halfbond or not np.array_equal(bonds.colors, colors):
        bonds.colors = colors
        # Refresh either native cylinders or the existing spacing renderer in
        # this frame. No new pseudobonds, labels or model IDs are created.
        group._update_graphics(group._ALL_CHANGE)
        group.redraw_needed(shape_changed=False)


def _update_distance_colors(session):
    from chimerax.atomic import PseudobondGroup
    for group in session.models.list(type=PseudobondGroup):
        _update_group_distance_colors(group)


def _ensure_distance_color_updates(session):
    # One session handler, regardless of the number of groups or replacements.
    # Rechecking at graphics updates also covers motion of an ancestor model.
    if getattr(session, '_xlms_distance_color_handler', None) is None:
        session._xlms_distance_color_handler = session.triggers.add_handler(
            'graphics update', lambda *_: _update_distance_colors(session))


def _install_distance_coloring(group, base_rgba, cutoff, violation_rgba, enabled):
    attributes = (('xlms_color_by_distance', bool(enabled), bool),
                  ('xlms_distance_cutoff', float(cutoff), float),
                  ('xlms_base_rgba', bytes(base_rgba).hex(), str),
                  ('xlms_violation_rgba', bytes(violation_rgba).hex(), str))
    for name, value, attr_type in attributes:
        group.register_attr(group.session, name, 'XLMS crosslinks', attr_type=attr_type)
        setattr(group, name, value)
    _ensure_distance_color_updates(group.session)
    _update_group_distance_colors(group)


def _style_labels(session, group, show_distances, label_height, label_size):
    from chimerax.core.commands import run
    from chimerax.label.label3d import labels_model
    distances = labels_model(group, create=True)
    options = []
    if label_height is not None:
        options.append(f'height {label_height:g}')
    if label_size is not None:
        options.append(f'size {label_size}')
    if options:
        run(session, f'label #{group.id_string} pseudobonds ' + ' '.join(options), log=False)
    distances.name = 'distances'
    distances.display = show_distances
    return distances


def _is_replaceable_group(model):
    """Global pseudobond groups, including legacy containers of such groups."""
    from chimerax.atomic import PseudobondGroup
    from chimerax.core.models import Model
    if isinstance(model, PseudobondGroup):
        return model.structure is None
    if type(model) is Model:
        children = model.child_models()
        return bool(children) and all(_is_replaceable_group(child) for child in children)
    return False


def _validate_output_model(session, output_model, replacing):
    """Resolve an exact destination without deleting or moving any models."""
    if output_model is None:
        return None
    text = str(output_model).strip()
    if not re.fullmatch(r'#?[1-9][0-9]*(?:\.[1-9][0-9]*)*', text):
        raise CrosslinkInputError("output_model must be an exact model ID such as '#2.1'.")
    destination = tuple(int(n) for n in text.lstrip('#').split('.'))
    models = {m.id_string: m for m in session.models.list()}
    target = '.'.join(str(n) for n in destination)
    parent_id = '.'.join(str(n) for n in destination[:-1])
    if parent_id and parent_id not in models:
        raise CrosslinkInputError(f'Parent model #{parent_id} must already exist for output_model=#{target}.')

    occupant = models.get(target)
    if occupant is not None and occupant not in replacing and _is_replaceable_group(occupant):
        replacing.append(occupant)

    def will_close(model_id):
        return any(model_id == m.id_string or model_id.startswith(m.id_string + '.')
                   for m in replacing)

    if parent_id and will_close(parent_id):
        raise CrosslinkInputError('The requested parent belongs to the crosslink group being replaced. '
                                  'Choose a parent outside that group.')
    if target in models and not will_close(target):
        raise CrosslinkInputError(f'output_model=#{target} is occupied by {models[target].name!r}. '
                                  'Only global pseudobond groups or their legacy containers can be replaced.')
    return destination


def _add(session, first, second, group_name, data_file, crosslink_color,
         show_distances, replace, radius, dashes, intramolecular=False,
         label_height=None, label_size=None, dash_spacing=None, output_model=None,
                         label_decimals=None, color_by_distance=True,
                         xl_cutoff=27.0, violation_color='red'):
    from chimerax.core.commands import ColorArg, StringArg, run
    from chimerax.core.errors import UserError

    try:
        if label_decimals is not None and (type(label_decimals) is not int or label_decimals < 0):
            raise CrosslinkInputError('label_decimals must be a nonnegative integer.')
        if not group_name.strip():
            raise CrosslinkInputError('group_name cannot be empty.')
        if not math.isfinite(radius) or radius <= 0 or not isinstance(dashes, int) or dashes < 0:
            raise CrosslinkInputError('radius must be positive and finite; dashes must be a nonnegative integer.')
        if label_height is not None and (not math.isfinite(label_height) or label_height <= 0):
            raise CrosslinkInputError('label_height must be positive and finite.')
        if label_size is not None and (not isinstance(label_size, int) or label_size <= 0):
            raise CrosslinkInputError('label_size must be a positive integer.')
        if dash_spacing is not None and (not math.isfinite(dash_spacing) or dash_spacing <= 0):
            raise CrosslinkInputError('dash_spacing must be positive and finite.')
        if not math.isfinite(xl_cutoff) or xl_cutoff <= 0:
            raise CrosslinkInputError('xl_cutoff must be positive and finite.')
        color_by_distance = parse_bool(color_by_distance)
        show_distances = parse_bool(show_distances)
        replace = parse_bool(replace)
        # Validate both colors before creating or replacing a model.
        base_rgba = _parse_color_rgba(session, crosslink_color, 'crosslink_color')
        violation_rgba = _parse_color_rgba(session, violation_color, 'violation_color')
        columns, pairs = read_pairs(data_file, intramolecular)
        chain1 = _chain_residues(session, first)
        chain2 = chain1 if intramolecular else _chain_residues(session, second)
        if not intramolecular and chain1[1] == chain2[1]:
            raise CrosslinkInputError('Intermolecular endpoints must specify different chains or models.')
        prepared, skipped, seen, lengths = [], [], set(), []
        duplicates = 0
        for line, pos1, pos2 in pairs:
            atom1, reason1 = _ca_atom(chain1, pos1)
            atom2, reason2 = _ca_atom(chain2, pos2)
            if reason1 or reason2:
                skipped.append(f'line {line}: ' + '; '.join(r for r in (reason1, reason2) if r))
                continue
            if atom1 is atom2:
                skipped.append(f'line {line}: both endpoints resolve to the same atom')
                continue
            key = tuple(sorted((atom1.atomspec, atom2.atomspec)))
            if key in seen:
                duplicates += 1
                continue
            seen.add(key)
            prepared.append((atom1.atomspec, atom2.atomspec))
            lengths.append(math.dist(atom1.scene_coord, atom2.scene_coord))
        if not prepared:
            details = '; '.join(skipped[:3])
            raise CrosslinkInputError(
                f'No drawable crosslinks: {len(skipped)} unresolved rows. '
                + details)
        if dash_spacing is not None:
            for length in lengths:
                dashes_for_spacing(length, dash_spacing)
    except (CrosslinkInputError, ValueError) as error:
        raise UserError(str(error)) from error

    # Explicit IDs identify the replacement target, independent of its name.
    # Without an output_model, preserve the case-insensitive name lookup.
    named = [m for m in session.models.list()
             if m.name.casefold() == group_name.casefold() and _is_replaceable_group(m)]
    existing = list(named) if output_model is None else []
    try:
        destination = _validate_output_model(session, output_model, existing)
    except CrosslinkInputError as error:
        raise UserError(str(error)) from error
    if existing and not replace:
        raise UserError('A pseudobond group already exists at the requested ID or name. '
                        'Choose another target or add replace=true.')
    # Global pseudobond names are unique. Never append to, or remove, a group
    # at a different ID merely because it has the newly requested name.
    if output_model is not None:
        conflicts = [m for m in named if not any(
            m is old or m.id_string.startswith(old.id_string + '.') for old in existing)]
        if conflicts:
            raise UserError(f'group_name={group_name!r} is already used at #{conflicts[0].id_string}. '
                            'Choose a different group_name or use that output_model.')
    if existing:
        session.models.close(existing)
    quoted_name = StringArg.unparse(group_name)
    quoted_color = StringArg.unparse(crosslink_color)
    group = None
    try:
        for spec1, spec2 in prepared:
            # One real pseudobond group owns every link and distance monitor.
            count = 0 if dash_spacing is not None else dashes
            command = (f'pbond {spec1} {spec2} name {quoted_name} global '
                       f'color {quoted_color} showDist true radius {radius:g} dashes {count}')
            run(session, command, log=False)
            if group is None:
                group = session.pb_manager.get_group(group_name, create=False)
                if group is None:
                    raise UserError('pbond did not create the requested global pseudobond group.')
        if destination is not None:
            # assign_id reparents the group and updates every child label ID.
            session.models.assign_id(group, destination)
        label_models = [_style_labels(session, group, show_distances, label_height, label_size)]
        _install_distance_coloring(group, base_rgba, xl_cutoff,
                                   violation_rgba, color_by_distance)
        if dash_spacing is not None:
            _install_spacing_renderer(group, dash_spacing)
        if label_decimals is not None:
            # Native distance-monitor formatting persists as distances update.
            # ChimeraX applies this preference to all monitored distances.
            session.pb_dist_monitor.decimal_places = label_decimals
    except Exception:
        if group is not None:
            session.models.close([group])
        failed = session.pb_manager.get_group(group_name, create=False)
        if failed is not None:
            session.models.close([failed])
        raise
    label_ids = ', '.join('#' + model.id_string for model in label_models)
    session.logger.info(
        f'{group_name}: created {len(prepared)} CA-to-CA crosslinks from {columns[0]} / {columns[1]}; '
        f'skipped {len(skipped)} unresolved rows, {duplicates} duplicate pairs. '
        f'Distance labels {"shown" if show_distances else "hidden"} in {label_ids}.'
    )
    session.logger.info(
        f'Initial CA-to-CA distances (Angstroms): min={min(lengths):.2f}, '
        f'mean={sum(lengths) / len(lengths):.2f}, max={max(lengths):.2f}. '
        f'Crosslink model #{group.id_string}.'
    )
    if dash_spacing is not None:
        session.logger.info(f'Target dash-plus-gap repeat spacing: {dash_spacing:g} Angstroms '
                            '(approximate; one crosslink model with one distances child).')
    if color_by_distance:
        session.logger.info(f'Live crosslink coloring: distances > {xl_cutoff:g} Angstroms use {violation_color}.')
    if skipped:
        session.logger.warning('Unresolved crosslinks (up to 20):\n' + '\n'.join(skipped[:20]))
    return group


def add_crosslinks_inter(session, protein_1, protein_2, group_name, data_file,
                         crosslink_color='blue', show_distances=False,
                         replace=True, radius=0.075, dashes=6,
                         label_height=None, label_size=None, dash_spacing=None, output_model=None,
                         label_decimals=None, color_by_distance=True,
                         xl_cutoff=27.0, violation_color='red'):
    return _add(session, protein_1, protein_2, group_name, data_file,
                crosslink_color, show_distances, replace, radius, dashes,
                label_height=label_height, label_size=label_size, dash_spacing=dash_spacing, output_model=output_model,
                label_decimals=label_decimals, color_by_distance=color_by_distance,
                xl_cutoff=xl_cutoff, violation_color=violation_color)


def add_crosslinks_intra(session, protein, group_name, data_file,
                                 crosslink_color='blue', show_distances=False,
                                 replace=True,
                                 radius=0.075, dashes=6,
                                 label_height=None, label_size=None, dash_spacing=None, output_model=None,
                         label_decimals=None, color_by_distance=True,
                         xl_cutoff=27.0, violation_color='red'):
    return _add(session, protein, None, group_name, data_file,
                crosslink_color, show_distances, replace, radius, dashes,
                intramolecular=True, label_height=label_height, label_size=label_size, dash_spacing=dash_spacing, output_model=output_model,
                label_decimals=label_decimals, color_by_distance=color_by_distance,
                xl_cutoff=xl_cutoff, violation_color=violation_color)


# Keep the longer Python name available as well.
add_crosslinks_intramolecular = add_crosslinks_intra


def _inter_command(session, options):
    from chimerax.core.errors import UserError
    try:
        values = parse_options(options)
    except CrosslinkInputError as error:
        raise UserError(str(error)) from error
    return add_crosslinks_inter(session, **values)


def _intra_command(session, options):
    from chimerax.core.errors import UserError
    try:
        values = parse_options(options, intramolecular=True)
    except CrosslinkInputError as error:
        raise UserError(str(error)) from error
    return add_crosslinks_intra(session, **values)


def register_commands(session):
    from chimerax.core.commands import CmdDesc, RestOfLine, register
    for name, function in (('add_crosslinks_inter', _inter_command),
                           ('add_crosslinks_intra', _intra_command),
                           ('add_crosslinks_intramolecular', _intra_command)):
        description = CmdDesc(required=[('options', RestOfLine)],
                              synopsis='Read an XL-MS CSV and add a named CA pseudobond group')
        register(name, description, function, logger=session.logger)
    _restore_spacing_renderers(session)
    if not getattr(session, '_xlms_restore_handler', None):
        session._xlms_restore_handler = session.triggers.add_handler(
            'end restore session', lambda *_: _restore_spacing_renderers(session))
    session.logger.info('Registered add_crosslinks_inter and add_crosslinks_intra (alias: add_crosslinks_intramolecular).')


if 'session' in globals():
    register_commands(session)
