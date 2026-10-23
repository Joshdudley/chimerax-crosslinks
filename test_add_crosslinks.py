"""Local regression checks with a lightweight ChimeraX API test double.

These verify parsing, chain routing, generated commands and label visibility;
they do not replace an in-application ChimeraX graphics test.
"""
import csv
import importlib.util
import shlex
import sys
import tempfile
import types
import unittest
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('xlms_commands', HERE / 'add_crosslinks.py')
xl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(xl)


class UserError(ValueError):
    pass


class FakeStructure:
    def __init__(self, model, chain, positions, missing_ca=()):
        self.id_string = model
        self.name = 'E6AP PDB'
        self.children = []
        self.residues = [FakeResidue(model, chain, n, n not in missing_ca) for n in positions]


class FakeResidue:
    def __init__(self, model, chain, number, has_ca):
        self.chain_id, self.number = chain, number
        self.atom = types.SimpleNamespace(atomspec=f'#{model}/{chain}:{number}@CA',
                                          scene_coord=(float(number), 0., 0.)) if has_ca else None

    def find_atom(self, name):
        assert name == 'CA'
        return self.atom


class FakeGroup:
    def __init__(self, name):
        self.name, self.labels, self.commands = name, None, []
        self.id_string = '9'
        self.structure = None
        self.deleted = False
        self.drawings = []
        self.session = None
        self._ALL_CHANGE = 255
        self._cylinder_sides = 12
        self._pbond_drawing = None
        self.pseudobonds = FakeBonds([])

    def register_attr(self, *args, **kwargs):
        pass

    def _update_graphics(self, changes=255):
        if self._pbond_drawing is None:
            self._pbond_drawing = FakeDrawing('native', None, None)
            self.add_drawing(self._pbond_drawing)
        self._pbond_drawing.visible_bonds = self.pseudobonds

    def _pseudobonds_planes_pick(self, planes):
        return list(self._pbond_drawing.visible_bonds.rows)

    def _update_positions(self, bonds, atoms):
        return [(tuple(a), tuple(b)) for a, b in zip(atoms[0].pb_scene_coords, atoms[1].pb_scene_coords)]

    def add_drawing(self, drawing):
        self.drawings.append(drawing)

    def remove_drawing(self, drawing):
        self.drawings.remove(drawing)

    def redraw_needed(self, **kwargs):
        pass


class FakeDrawing:
    def __init__(self, name, pick, picks):
        self.name, self.display = name, True

    def set_geometry(self, *geometry):
        self.geometry = geometry


class FakeBonds:
    def __init__(self, rows, colors=None):
        self.rows = rows
        self._colors = None if colors is None else np.array(colors, dtype=np.uint8)
        self.halfbonds = False
        self.color_writes = 0

    def __len__(self):
        return len(self.rows)

    @property
    def colors(self):
        if self._colors is None or self._colors.shape != (len(self), 4):
            self._colors = np.tile([0, 0, 255, 255], (len(self), 1)).astype(np.uint8)
        return self._colors.copy()

    @colors.setter
    def colors(self, value):
        self._colors = np.broadcast_to(np.asarray(value, dtype=np.uint8), (len(self), 4)).copy()
        self.color_writes += 1

    @property
    def lengths(self):
        first, second = self.atoms
        return np.linalg.norm(first.scene_coords - second.scene_coords, axis=1)

    @property
    def atoms(self):
        coords = [np.array([row[i] for row in self.rows], dtype=float).reshape((-1, 3)) for i in (0, 1)]
        return tuple(types.SimpleNamespace(pb_scene_coords=c, scene_coords=c) for c in coords)

    def filter(self, mask):
        return FakeBonds([row for row, keep in zip(self.rows, mask) if keep], self.colors[mask])

    @property
    def half_colors(self):
        return np.concatenate((self.colors, self.colors), axis=0)


class FakeContainer:
    def __init__(self, name, session):
        self.name, self.session, self.children = name, session, []
        self.deleted = False

    def child_models(self):
        return self.children


class FakeTriggers:
    def __init__(self):
        self.handlers = {}

    def add_handler(self, event, callback):
        self.handlers.setdefault(event, []).append(callback)
        return callback

    def fire(self, event):
        for callback in tuple(self.handlers.get(event, ())):
            callback(event, None)


class FakeSession:
    def __init__(self, *structures):
        self.groups, self.commands, self.messages = {}, [], []
        self.triggers = FakeTriggers()
        self.containers = []
        self.next_id = 9
        self.structures = structures
        self.models = types.SimpleNamespace(list=self.list_models, close=self.close, add=self.add, assign_id=self.assign_id)
        self.pb_manager = types.SimpleNamespace(get_group=self.get_group)
        self.logger = types.SimpleNamespace(info=self.messages.append, warning=self.messages.append)
        self.pb_dist_monitor = types.SimpleNamespace(decimal_places=3)
        self.fail_at = None

    def get_group(self, name, create=True):
        if name not in self.groups and create:
            group = self.groups[name] = FakeGroup(name)
            group.session = self
            group.id_string = str(self.next_id)
            self.next_id += 1
        return self.groups.get(name)

    def list_models(self, type=None):
        models = list(self.structures) + self.containers + list(self.groups.values())
        return models if type is None else [m for m in models if isinstance(m, type)]

    def add(self, models, parent=None):
        for model in models:
            if isinstance(model, FakeContainer) and model not in self.containers:
                self.containers.append(model)
                model.id_string = str(self.next_id)
                self.next_id += 1
            if parent is not None:
                parent.children.append(model)
                model.id_string = f'{parent.id_string}.{len(parent.children)}'

    def assign_id(self, model, destination):
        target = '.'.join(str(n) for n in destination)
        if any(m.id_string == target and m is not model for m in self.list_models()):
            raise ValueError('occupied model ID')
        parent_id = '.'.join(str(n) for n in destination[:-1])
        parent = next((m for m in self.list_models() if m.id_string == parent_id), None)
        old_parent = getattr(model, 'parent', None)
        if old_parent is not None:
            old_parent.children.remove(model)
        model.id_string = target
        model.parent = parent
        if parent is not None:
            parent.children.append(model)
        if model.labels is not None:
            model.labels.id_string = target + '.1'

    def close(self, groups):
        for group in groups:
            if group.deleted:
                continue
            if isinstance(group, FakeContainer):
                self.close(group.children)
                self.containers.remove(group)
            else:
                self.groups.pop(group.name, None)
            parent = getattr(group, 'parent', None)
            if parent is not None:
                parent.children.remove(group)
            group.deleted = True


def fake_run(session, command, log=False):
    if session.fail_at == len(session.commands):
        raise RuntimeError('simulated pbond failure')
    tokens = shlex.split(command)
    if tokens[0] == 'label':
        group = next(g for g in session.groups.values() if '#' + g.id_string == tokens[1])
        group.label_command = command
        return
    assert tokens[0] == 'pbond'
    assert tokens[1].endswith('@CA') and tokens[2].endswith('@CA')
    assert 'global' in tokens
    assert tokens[tokens.index('showDist') + 1] == 'true'
    name = tokens[tokens.index('name') + 1]
    group = session.pb_manager.get_group(name)
    if group.labels is None:
        group.labels = types.SimpleNamespace(name='labels', display=True, id_string=group.id_string + '.1')
    atoms = {r.atom.atomspec: r.atom for structure in session.structures for r in structure.residues if r.atom}
    group.pseudobonds.rows.append((atoms[tokens[1]].scene_coord, atoms[tokens[2]].scene_coord))
    colors = group.pseudobonds.colors
    colors[-1] = fake_color_parse(tokens[tokens.index('color') + 1], session)[0].uint8x4()
    group.pseudobonds.colors = colors
    group.commands.append(command)
    session.commands.append(command)


def fake_color_parse(text, session):
    if 'invalid-color' in text:
        raise ValueError('invalid color')
    values = {'blue': (0, 0, 255, 255), 'red': (255, 0, 0, 255),
              'yellow': (255, 255, 0, 255), 'green': (0, 255, 0, 255)}
    rgba = values.get(shlex.split(text)[0], (0, 0, 255, 255))
    return types.SimpleNamespace(uint8x4=lambda: rgba), text, ''


for name in ('chimerax', 'chimerax.core', 'chimerax.core.commands',
             'chimerax.core.errors', 'chimerax.core.models', 'chimerax.atomic', 'chimerax.atomic.structure', 'chimerax.surface', 'chimerax.label', 'chimerax.label.label3d'):
    sys.modules[name] = types.ModuleType(name)
cmd = sys.modules['chimerax.core.commands']
cmd.run = fake_run
cmd.ColorArg = types.SimpleNamespace(parse=fake_color_parse)
cmd.StringArg = types.SimpleNamespace(unparse=shlex.quote)
sys.modules['chimerax.core.errors'].UserError = UserError
sys.modules['chimerax.core.models'].Model = FakeContainer
sys.modules['chimerax.atomic'].AtomicStructure = FakeStructure
sys.modules['chimerax.atomic'].PseudobondGroup = FakeGroup
graphics = sys.modules['chimerax.atomic.structure']
graphics.BondsDrawing = FakeDrawing
graphics.PickedPseudobond = graphics.PickedPseudobonds = object
graphics._selected_bond_cylinders = lambda bonds: [False] * (len(bonds.rows) * 2)
sys.modules['chimerax.surface'].dashed_cylinder_geometry = lambda n, height, nc: (n, height, nc)
sys.modules['chimerax.label.label3d'].labels_model = lambda group, create=False: group.labels


class CrosslinkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def fixture(self, text):
        p = Path(self.temp.name) / 'pairs.csv'
        p.write_text(text)
        return str(p)

    def test_quoted_options_and_protein_spec(self):
        values = xl.parse_options("protein_1='E6AP, model #1.1/A' protein_2='CaMKII_WT, model #1.2/B' "
                                  "group_name='my links' data_file='/a path/input.csv' show_distances=false")
        self.assertIs(values['show_distances'], False)
        self.assertEqual(xl.parse_protein(values['protein_1']), ('E6AP', '1.1', 'A'))
        self.assertEqual(values['data_file'], '/a path/input.csv')
        with self.assertRaises(xl.CrosslinkInputError):
            xl.parse_options('protein=x group_name=x data_file=x bogus=true', True)

    def test_inter_routes_by_columns_not_protein_label(self):
        path = self.fixture('E6AP,CaMKII_D135N,self_dimer\nK774,K199,false\n')
        s = FakeSession(FakeStructure('1.1', 'A', [774]), FakeStructure('1.2', 'B', [199]))
        group = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII_WT, model #1.2/B', 'links', path)
        self.assertIn('pbond #1.1/A:774@CA #1.2/B:199@CA', s.commands[0])
        self.assertEqual(group.labels.name, 'distances')
        self.assertFalse(group.labels.display)

    def test_intra_draws_flagged_rows_and_shows_labels(self):
        path = self.fixture('Residue1,Residue2,self_dimer\nK27,E33,true\nE33,K27,false\nK27,K88,false\n')
        s = FakeSession(FakeStructure('1.1', 'A', [27, 33]))
        group = xl.add_crosslinks_intra(s, 'E6AP, model #1.1/A', 'intra', path, show_distances=True)
        self.assertEqual(len(s.commands), 1)
        self.assertTrue(group.labels.display)
        self.assertIn('pbond #1.1/A:27@CA #1.1/A:33@CA', s.commands[0])
        self.assertTrue(any('1 unresolved rows, 1 duplicate pairs' in m for m in s.messages))

    def test_same_protein_different_copies(self):
        path = self.fixture('Residue1,Residue2,self_dimer\nK27,E33,true\n')
        s = FakeSession(FakeStructure('1.1', 'A', [27]), FakeStructure('1.3', 'C', [33]))
        xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'E6AP, model #1.3/C', 'self', path)
        self.assertIn('pbond #1.1/A:27@CA #1.3/C:33@CA', s.commands[0])

    def test_missing_extension_and_legacy_two_columns(self):
        path = self.fixture('Residue1,Residue2\nK10,20\n')
        columns, rows = xl.read_pairs(path[:-4], True)
        self.assertEqual(columns, ('Residue1', 'Residue2'))
        self.assertEqual(rows, [(2, 10, 20)])

    def test_validation_does_not_replace_good_group(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10]), FakeStructure('1.2', 'B', [20]))
        original = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path)
        with self.assertRaises(UserError):
            xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path, replace=False)
        with self.assertRaises(UserError):
            xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path,
                                    replace=True, crosslink_color='invalid-color')
        self.assertIs(s.groups['links'], original)
        replacement = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path,
                                             show_distances=True)
        self.assertIsNot(replacement, original)
        self.assertEqual(len(replacement.commands), 1)
        self.assertIsNot(replacement.labels, original.labels)
        self.assertTrue(replacement.labels.display)

    def test_intra_replaces_group_by_default(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10, 20]))
        original = xl.add_crosslinks_intra(s, 'E6AP, model #1.1/A', 'intra', path)
        replacement = xl.add_crosslinks_intra(s, 'E6AP, model #1.1/A', 'intra', path)
        self.assertIsNot(replacement, original)
        self.assertEqual(len(s.groups), 1)
        self.assertEqual(len(replacement.commands), 1)

    def test_partial_group_cleaned_on_command_failure(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\nK11,K20\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10, 11]), FakeStructure('1.2', 'B', [20]))
        s.fail_at = 1
        with self.assertRaises(RuntimeError):
            xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path)
        self.assertNotIn('links', s.groups)

    def test_optional_label_controls_and_distance_summary(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\nK12,K20\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10, 12]), FakeStructure('1.2', 'B', [20]))
        options = xl.parse_options(f"protein_1='E6AP, model #1.1/A' protein_2='CaMKII, model #1.2/B' "
                                   f"group_name=links data_file='{path}' label_height=1.5 label_size=64")
        group = xl.add_crosslinks_inter(s, **options)
        self.assertEqual(group.label_command, 'label #9 pseudobonds height 1.5 size 64')
        self.assertFalse(group.labels.display)
        self.assertTrue(any('min=8.00, mean=9.00, max=10.00' in m for m in s.messages))
        with self.assertRaises(xl.CrosslinkInputError):
            xl.parse_options('protein=x group_name=x data_file=x label_size=0', True)

    def test_example_csv_files(self):
        paths = sorted((HERE / 'examples').glob('*.csv'))
        self.assertEqual(len(paths), 2)
        counts = []
        for path in paths:
            intra = 'intramolecular' in path.name
            columns, pairs = xl.read_pairs(str(path), intra)
            counts.append(len(pairs))
            if intra:
                positions = sorted({n for _, a, b in pairs for n in (a, b)})
                s = FakeSession(FakeStructure('1.1', 'A', positions))
                xl.add_crosslinks_intra(s, 'protein, model #1.1/A', path.stem, str(path))
                expected = {tuple(sorted((a, b))) for _, a, b in pairs if a != b}
            else:
                s = FakeSession(FakeStructure('1.1', 'A', sorted({a for _, a, b in pairs})),
                                FakeStructure('1.2', 'B', sorted({b for _, a, b in pairs})))
                xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', path.stem, str(path))
                expected = {(a, b) for _, a, b in pairs}
            self.assertEqual(len(s.commands), len(expected))
            print(path.name, ':', len(pairs), 'rows ->', len(s.commands), 'drawable unique pairs in test model')
        self.assertEqual(counts, [2, 2])

    def test_spacing_scales_counts_and_replaces_entire_group(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\nK10,K30\nK12,K22\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10, 12]), FakeStructure('1.2', 'B', [20, 22, 30]))
        root = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path,
                                      dash_spacing=1, label_height=1.5, show_distances=False)
        self.assertEqual(len(s.groups), 1)
        self.assertFalse(s.containers)
        self.assertEqual(set(root._xlms_renderer.drawings), {10, 20})
        self.assertEqual(len(root.commands), 3)
        self.assertFalse(root._pbond_drawing.display)
        self.assertEqual(root.labels.name, 'distances')
        self.assertFalse(root.labels.display)
        self.assertEqual(root._xlms_renderer.drawings[10].geometry, (5, 0.5, 12))
        self.assertEqual(len(root._xlms_renderer.drawings[10].visible_bonds.rows), 2)
        self.assertEqual(len(root._pseudobonds_planes_pick(None)), 3)
        replacement = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path,
                                             dash_spacing=0.5, show_distances=True)
        self.assertTrue(root.deleted)
        self.assertEqual(len(s.groups), 1)
        self.assertEqual(set(replacement._xlms_renderer.drawings), {20, 40})
        self.assertTrue(replacement.labels.display)
        # Switching back to fixed-count mode must remove the spacing-mode root.
        plain = xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path)
        self.assertTrue(replacement.deleted)
        self.assertEqual(len(s.groups), 1)
        self.assertEqual(plain.name, 'links')

    def test_spacing_failure_cleans_children(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\nK10,K30\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10]), FakeStructure('1.2', 'B', [20, 30]))
        s.fail_at = 1
        with self.assertRaises(RuntimeError):
            xl.add_crosslinks_inter(s, 'E6AP, model #1.1/A', 'CaMKII, model #1.2/B', 'links', path,
                                    dash_spacing=1)
        self.assertFalse(s.groups)
        self.assertFalse(s.containers)

    def test_spacing_renderer_updates_after_motion_and_hiding(self):
        group = FakeGroup('links')
        group.pseudobonds = FakeBonds([((0, 0, 0), (10, 0, 0))])
        xl._install_spacing_renderer(group, 1)
        self.assertEqual(set(group._xlms_renderer.drawings), {10})
        group.pseudobonds.rows[0] = ((0, 0, 0), (30, 0, 0))
        group._update_graphics()
        self.assertEqual(set(group._xlms_renderer.drawings), {30})
        self.assertEqual(group._xlms_renderer.drawings[30].positions, [((0., 0., 0.), (30., 0., 0.))])
        group.pseudobonds = FakeBonds([])
        group._update_graphics()
        self.assertFalse(group._xlms_renderer.drawings)
        self.assertEqual(len(group.drawings), 1)  # Hidden native drawing only.

    def test_saved_spacing_is_reinstalled_without_extra_models(self):
        s = FakeSession()
        group = s.get_group('links')
        group.xlms_dash_spacing = 1.0  # The attribute restored by ChimeraX.
        group.pseudobonds = FakeBonds([((0, 0, 0), (12, 0, 0))])
        xl._restore_spacing_renderers(s)
        renderer = group._xlms_renderer
        self.assertEqual(set(renderer.drawings), {12})
        xl._restore_spacing_renderers(s)
        self.assertIs(group._xlms_renderer, renderer)
        self.assertEqual(len(group.drawings), 2)  # One hidden native, one custom.
        self.assertEqual(len(s.groups), 1)

    def test_legacy_subgroups_are_replaced(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('1.1', 'A', [10, 20]))
        old = FakeContainer('links', s)
        s.add([old])
        child = s.get_group('links [dashes=10, owner=#9]')
        s.add([child], parent=old)
        new = xl.add_crosslinks_intra(s, 'E6AP, model #1.1/A', 'links', path, dash_spacing=1)
        self.assertTrue(old.deleted)
        self.assertTrue(child.deleted)
        self.assertEqual(list(s.groups), ['links'])
        self.assertFalse(s.containers)
        self.assertEqual(new.labels.name, 'distances')

    def test_output_model_nests_links_and_labels_and_replaces(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        values = xl.parse_options(f"protein='E6AP, model #2/A' group_name=Intra_E6AP "
                                  f"data_file='{path}' output_model='#2.1' dash_spacing=1 label_height=1.5", True)
        group = xl.add_crosslinks_intra(s, **values)
        self.assertEqual(group.id_string, '2.1')
        self.assertIs(group.parent, s.structures[0])
        self.assertEqual(group.labels.id_string, '2.1.1')
        self.assertIn('label #2.1 pseudobonds', group.label_command)
        self.assertEqual(set(group._xlms_renderer.drawings), {10})
        replacement = xl.add_crosslinks_intra(s, **values)
        self.assertTrue(group.deleted)
        self.assertEqual(replacement.id_string, '2.1')
        self.assertEqual(s.structures[0].children, [replacement])

    def test_inter_output_model_and_top_level_output(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10]), FakeStructure('3', 'B', [20]))
        group = xl.add_crosslinks_inter(s, 'E6AP, model #2/A', 'CaMKII, model #3/B', 'links', path,
                                        output_model='#2.4', show_distances=True)
        self.assertEqual(group.id_string, '2.4')
        self.assertEqual(group.labels.id_string, '2.4.1')
        group = xl.add_crosslinks_inter(s, 'E6AP, model #2/A', 'CaMKII, model #3/B', 'other_links', path,
                                        output_model='#8')
        self.assertEqual(group.id_string, '8')
        self.assertEqual(group.labels.id_string, '8.1')
        self.assertIsNone(group.parent)

    def test_output_model_errors_preserve_existing_links(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        original = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path, output_model='#2.1')
        other = s.get_group('unrelated')
        s.assign_id(other, (2, 2))
        # Parent must survive replacement; occupied IDs and malformed selectors
        # must be rejected before closing or creating any models.
        for value in ('#2', '#2.2', '#5.1', '#2.1.3', '#2.1/A', '#0', '#2.0', '#2,#3'):
            with self.subTest(value=value), self.assertRaises(UserError):
                xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path, output_model=value)
            self.assertIs(s.groups['links'], original)
            self.assertFalse(original.deleted)
            self.assertFalse(other.deleted)
        with self.assertRaises(UserError):
            xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path,
                                    output_model='#2.1', replace=False)
        self.assertFalse(original.deleted)

    def test_inter_case_change_replaces_at_requested_id(self):
        path = self.fixture('E6AP,CaMKII\nK10,K20\n')
        s = FakeSession(FakeStructure('4', 'A', [10]), FakeStructure('5', 'B', [20]))
        old_name = 'Inter_E6AP_to_CaMKII_Group_1'
        new_name = 'Inter_E6AP_to_CaMKII_group_1'
        old = xl.add_crosslinks_inter(s, 'E6AP, model #4/A', 'CaMKII_WT, model #5/B', old_name, path,
                                      output_model='#5.1', dash_spacing=1)
        # Bad input and replace=false must preserve the existing case variant.
        with self.assertRaises(UserError):
            xl.add_crosslinks_inter(s, 'E6AP, model #4/A', 'CaMKII_WT, model #5/B', new_name, path,
                                    output_model='#5.1', replace=False)
        with self.assertRaises(UserError):
            xl.add_crosslinks_inter(s, 'E6AP, model #4/A', 'CaMKII_WT, model #5/B', new_name,
                                    path + '.missing', output_model='#5.1')
        self.assertFalse(old.deleted)
        new = xl.add_crosslinks_inter(s, 'E6AP, model #4/A', 'CaMKII_WT, model #5/B', new_name, path,
                                      output_model='#5.1', crosslink_color='yellow', show_distances=True,
                                      label_height=1.5, label_size=64, label_decimals=1,
                                      dash_spacing=1, radius=0.2)
        self.assertTrue(old.deleted)
        self.assertEqual(list(s.groups), [new_name])
        self.assertEqual(new.id_string, '5.1')
        self.assertEqual(new.labels.id_string, '5.1.1')
        self.assertIsNot(new.labels, old.labels)
        self.assertTrue(new.labels.display)
        self.assertEqual(s.structures[1].children, [new])

    def test_intra_case_change_replaces_without_output_id(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        old = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'Intra_E6AP', path)
        new = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'intra_e6ap', path)
        self.assertTrue(old.deleted)
        self.assertEqual(list(s.groups), ['intra_e6ap'])
        self.assertEqual(new.name, 'intra_e6ap')

    def test_case_match_does_not_replace_unrelated_target(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        old = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'Intra_E6AP', path)
        other = s.get_group('Unrelated')
        s.assign_id(other, (2, 1))
        with self.assertRaises(UserError):
            xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'intra_e6ap', path, output_model='#2.1')
        self.assertFalse(old.deleted)
        self.assertFalse(other.deleted)

    def test_output_id_replaces_completely_different_name(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        old = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'Old_links', path, output_model='#2.1')
        with self.assertRaises(UserError):
            xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'New_links', path,
                                    output_model='#2.1', replace=False)
        self.assertFalse(old.deleted)
        new = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'New_links', path, output_model='#2.1')
        self.assertTrue(old.deleted)
        self.assertEqual(list(s.groups), ['New_links'])
        self.assertEqual(new.id_string, '2.1')
        self.assertEqual(new.labels.id_string, '2.1.1')
        self.assertIsNot(new.labels, old.labels)

    def test_output_id_name_collision_is_checked_before_deleting(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        first = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'First', path, output_model='#2.1')
        second = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'Second', path, output_model='#2.2')
        with self.assertRaises(UserError):
            xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'Second', path, output_model='#2.1')
        self.assertFalse(first.deleted)
        self.assertFalse(second.deleted)
        self.assertEqual(len(s.groups), 2)

    def test_output_id_replaces_legacy_group_regardless_of_name(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        old = FakeContainer('Old container', s)
        s.add([old])
        child = s.get_group('Old count bucket')
        s.add([child], parent=old)
        old_id = '#' + old.id_string
        new = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'New_links', path, output_model=old_id)
        self.assertTrue(old.deleted)
        self.assertTrue(child.deleted)
        self.assertEqual('#' + new.id_string, old_id)

    def test_distance_precision_for_both_commands(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]), FakeStructure('3', 'B', [20]))
        options = xl.parse_options(f"protein='E6AP, model #2/A' group_name=intra data_file='{path}' "
                                   "label_decimals=1 output_model='#2.1'", True)
        group = xl.add_crosslinks_intra(s, **options)
        self.assertEqual(s.pb_dist_monitor.decimal_places, 1)
        self.assertEqual(group.labels.id_string, '2.1.1')
        self.assertFalse(group.labels.display)
        # Omitting precision preserves the native preference.
        xl.add_crosslinks_inter(s, 'E6AP, model #2/A', 'CaMKII, model #3/B', 'inter', path)
        self.assertEqual(s.pb_dist_monitor.decimal_places, 1)
        xl.add_crosslinks_inter(s, 'E6AP, model #2/A', 'CaMKII, model #3/B', 'inter', path,
                                label_decimals=0, show_distances=True)
        self.assertEqual(s.pb_dist_monitor.decimal_places, 0)
        # Python callers use the same option.
        xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'intra', path, label_decimals=2)
        self.assertEqual(s.pb_dist_monitor.decimal_places, 2)

    def test_bad_precision_preserves_existing_group(self):
        path = self.fixture('Residue1,Residue2\nK10,K20\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 20]))
        original = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'intra', path)
        for value in (-1, 1.5, True):
            with self.assertRaises(UserError):
                xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'intra', path, label_decimals=value)
            self.assertIs(s.groups['intra'], original)
        for extra in ('label_decimals=-1', 'label_decimals=1.5',
                      'label_decimals=1 label_decimals=2', 'label_sig_figs=1', 'label_decimal_places=1'):
            with self.assertRaises(xl.CrosslinkInputError):
                xl.parse_options('protein=x group_name=x data_file=x ' + extra, True)
        self.assertFalse(original.deleted)

    def test_live_colors_cross_threshold_in_both_directions(self):
        path = self.fixture('E6AP,CaMKII\nK10,K36\nK10,K37\nK10,K38\n')
        for spacing in (None, 1.0):
            with self.subTest(dash_spacing=spacing):
                s = FakeSession(FakeStructure('4', 'A', [10]), FakeStructure('5', 'B', [36, 37, 38]))
                group = xl.add_crosslinks_inter(s, 'E6AP, model #4/A', 'CaMKII, model #5/B',
                                                'links', path, crosslink_color='yellow',
                                                output_model='#5.1', dash_spacing=spacing, label_decimals=1)
                yellow, red = [255, 255, 0, 255], [255, 0, 0, 255]
                self.assertEqual(group.pseudobonds.colors.tolist(), [yellow, yellow, red])
                # Move the endpoint across the cutoff; no command is rerun.
                group.pseudobonds.rows[0] = ((10, 0, 0), (37.001, 0, 0))
                group.pseudobonds.rows[2] = ((10, 0, 0), (36, 0, 0))
                s.triggers.fire('graphics update')
                self.assertEqual(group.pseudobonds.colors.tolist(), [red, yellow, yellow])
                self.assertEqual(len(s.groups), 1)
                self.assertEqual(group.labels.id_string, '5.1.1')
                if spacing is not None:
                    # All have the same rounded dash count; their colors must
                    # still differ inside that one internal drawing.
                    drawing = group._xlms_renderer.drawings[28]
                    self.assertEqual(drawing.colors[:2].tolist(), [red, yellow])
                writes = group.pseudobonds.color_writes
                s.triggers.fire('graphics update')
                self.assertEqual(group.pseudobonds.color_writes, writes)

    def test_live_colors_for_intra_and_custom_settings(self):
        path = self.fixture('Residue1,Residue2\nK10,K38\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 38]))
        options = xl.parse_options(f"protein='E6AP, model #2/A' group_name=links data_file='{path}' "
                                   "xl_cutoff=25 violation_color=green color_by_distance=true", True)
        group = xl.add_crosslinks_intra(s, **options)
        self.assertEqual(group.pseudobonds.colors.tolist(), [[0, 255, 0, 255]])
        group.pseudobonds.rows[0] = ((10, 0, 0), (35, 0, 0))
        s.triggers.fire('graphics update')
        self.assertEqual(group.pseudobonds.colors.tolist(), [[0, 0, 255, 255]])
        replacement = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path,
                                              color_by_distance=False, crosslink_color='yellow')
        s.triggers.fire('graphics update')
        self.assertEqual(replacement.pseudobonds.colors.tolist(), [[255, 255, 0, 255]])
        self.assertEqual(len(s.triggers.handlers['graphics update']), 1)
        s.models.close([replacement])
        s.triggers.fire('graphics update')  # Deleted groups are no longer visited.

    def test_live_colors_restore_and_no_untagged_groups(self):
        s = FakeSession()
        unrelated = s.get_group('unrelated')
        unrelated.pseudobonds = FakeBonds([((0, 0, 0), (100, 0, 0))])
        restored = s.get_group('restored')
        restored.pseudobonds = FakeBonds([((0, 0, 0), (28, 0, 0))])
        restored.xlms_color_by_distance = True
        restored.xlms_distance_cutoff = 27.0
        restored.xlms_base_rgba = '0000ffff'
        restored.xlms_violation_rgba = 'ff0000ff'
        xl._restore_spacing_renderers(s)
        self.assertEqual(restored.pseudobonds.colors.tolist(), [[255, 0, 0, 255]])
        self.assertEqual(unrelated.pseudobonds.colors.tolist(), [[0, 0, 255, 255]])
        restored.pseudobonds.rows.clear()
        s.triggers.fire('graphics update')  # Empty groups are safe, too.

    def test_invalid_live_options_preserve_existing_group(self):
        path = self.fixture('Residue1,Residue2\nK10,K38\n')
        s = FakeSession(FakeStructure('2', 'A', [10, 38]))
        old = xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path)
        for option in ({'xl_cutoff': 0}, {'xl_cutoff': float('nan')},
                       {'violation_color': 'invalid-color'}, {'color_by_distance': 'maybe'}):
            with self.assertRaises(UserError):
                xl.add_crosslinks_intra(s, 'E6AP, model #2/A', 'links', path, **option)
            self.assertFalse(old.deleted)
        for value in ('0', '-1', 'nan', 'inf'):
            with self.assertRaises(xl.CrosslinkInputError):
                xl.parse_options('protein=x group_name=x data_file=x xl_cutoff=' + value, True)

    def test_spacing_rounding_and_limits(self):
        self.assertEqual(xl.dashes_for_spacing(30, 1), 30)
        self.assertEqual(xl.dashes_for_spacing(60, 1), 60)
        self.assertEqual(xl.dashes_for_spacing(0.1, 1), 2)
        self.assertEqual(xl.dashes_for_spacing(11, 1), 12)
        with self.assertRaises(xl.CrosslinkInputError):
            xl.dashes_for_spacing(100, 0.00001)
        with self.assertRaises(xl.CrosslinkInputError):
            xl.parse_options('protein=x group_name=x data_file=x dash_spacing=0', True)


if __name__ == '__main__':
    unittest.main()
