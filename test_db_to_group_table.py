#!/usr/bin/env python3
"""Checks of the density denominator's reject set in db_to_group_table.py.

Synthetic databases only. Run from the repository root with
    python3 -m unittest test_db_to_group_table
(standard library only; pytest also collects it).
"""
import csv
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest

BRIDGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'db_to_group_table.py')

# analysed_time as TurtleWave 4.3 and 4.4+ write it (dbwrite.store_analysed_time).
ANALYSED_43 = """CREATE TABLE analysed_time (
    subject TEXT NOT NULL, stage TEXT NOT NULL, reject_artifacts INTEGER NOT NULL,
    reject_arousals INTEGER NOT NULL, analysed_seconds REAL NOT NULL,
    artefact_seconds_excluded REAL, epoch_length REAL, source TEXT, annotation_file TEXT,
    turtlewave_version TEXT, processing_timestamp TEXT,
    PRIMARY KEY (subject, stage, reject_artifacts, reject_arousals))"""
ANALYSED_44 = """CREATE TABLE analysed_time (
    subject TEXT NOT NULL, stage TEXT NOT NULL, reject_types TEXT NOT NULL,
    reject_artifacts INTEGER NOT NULL, reject_arousals INTEGER NOT NULL,
    analysed_seconds REAL NOT NULL, artefact_seconds_excluded REAL, epoch_length REAL,
    source TEXT, annotation_file TEXT, turtlewave_version TEXT, processing_timestamp TEXT,
    PRIMARY KEY (subject, stage, reject_types))"""

# Two E1 events and one E2 event of one 12-15 Hz spindle run stored under NREM2NREM3.
EVENTS = [('E1', 1.0, 40.0), ('E1', 2.0, 60.0), ('E2', 1.5, 30.0)]


def make_subject(root, name, analysed):
    """analysed: {reject set: {stage: seconds}}; a set given as (reject_artifacts,
    reject_arousals) writes TurtleWave 4.3's table, as text the reject_types table of 4.4."""
    folder = os.path.join(root, name, 'wonambi')
    os.makedirs(folder)
    con = sqlite3.connect(os.path.join(folder, 'neural_events.db'))
    con.execute('CREATE TABLE events (channel TEXT, event_type TEXT, method TEXT, stage TEXT, '
                'freq_lower REAL, freq_upper REAL, duration REAL, peak2peak_amp REAL)')
    con.executemany("INSERT INTO events VALUES (?, 'spindle', 'm', 'NREM2NREM3', 12, 15, ?, ?)",
                    EVENTS)
    con.execute('CREATE TABLE processing_status (channel TEXT, event_type TEXT, method TEXT, '
                'freq_lower REAL, freq_upper REAL, stage TEXT, success INTEGER)')
    con.executemany("INSERT INTO processing_status VALUES (?, 'spindle', 'm', 12, 15, "
                    "'NREM2NREM3', 1)", [('E1',), ('E2',), ('E3',)])
    old = all(isinstance(k, tuple) for k in analysed)
    con.execute(ANALYSED_43 if old else ANALYSED_44)
    for key, seconds in analysed.items():
        for stage, sec in seconds.items():
            if old:
                con.execute('INSERT INTO analysed_time (subject, stage, reject_artifacts, '
                            "reject_arousals, analysed_seconds, artefact_seconds_excluded, "
                            "epoch_length, source, turtlewave_version) "
                            "VALUES (?, ?, ?, ?, ?, 0, 30, 'detection', '4.3.0')",
                            (name, stage, *key, sec))
            else:
                flags = (int('Artefact' in key), int('Arousal' in key))
                con.execute('INSERT INTO analysed_time (subject, stage, reject_types, '
                            'reject_artifacts, reject_arousals, analysed_seconds, '
                            'artefact_seconds_excluded, epoch_length, source, '
                            "turtlewave_version) VALUES (?, ?, ?, ?, ?, ?, 0, 30, "
                            "'detection', '4.6.2')", (name, stage, key, *flags, sec))
    con.commit()
    con.close()
    with open(os.path.join(folder, f'{name}_clean_rebuilt.xml'), 'w') as f:
        f.write('<annotations><rater><stages><epoch><epoch_start>0</epoch_start><epoch_end>30'
                '</epoch_end><stage>NREM2</stage></epoch></stages></rater></annotations>')


class RejectSet(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = os.path.join(self.tmp.name, 'in')
        self.out = os.path.join(self.tmp.name, 'out')

    def tearDown(self):
        self.tmp.cleanup()

    def bridge(self, *extra):
        return subprocess.run([sys.executable, BRIDGE, self.root, self.out,
                               '--bands', 'spindle_12-15', *extra],
                              capture_output=True, text=True)

    def ok(self, *extra):
        done = self.bridge(*extra)
        self.assertEqual(done.returncode, 0, done.stderr + done.stdout)
        with open(os.path.join(self.out, 'eventStat_density_spindle_12-15.csv')) as f:
            dens = {r['Subject']: r for r in csv.DictReader(f)}
        with open(os.path.join(self.out, 'recordings.csv')) as f:
            recs = {r['Subject']: r for r in csv.DictReader(f)}
        return dens, recs

    def refused(self, pattern, *extra):
        done = self.bridge(*extra)
        self.assertNotEqual(done.returncode, 0, done.stdout)
        self.assertRegex(done.stderr, pattern)
        self.assertFalse(os.path.exists(os.path.join(self.out, 'recordings.csv')))

    def test_43_artefact_only_table_divides_by_its_time(self):
        # What TurtleWave 4.3 detection stored: artefacts removed, arousals kept.
        make_subject(self.root, 's01', {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0}})
        dens, recs = self.ok()
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 40)
        self.assertAlmostEqual(float(dens['s01']['E3']), 0.0)
        self.assertEqual(recs['s01']['reject_types'], 'Artefact')
        self.assertAlmostEqual(float(recs['s01']['denominator_minutes']), 40.0)

    def test_43_both_flags_table_uses_both(self):
        make_subject(self.root, 's01', {(1, 1): {'NREM2': 900.0, 'NREM3': 300.0}})
        dens, recs = self.ok()
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 20)
        self.assertEqual(recs['s01']['reject_types'], 'Arousal,Artefact')

    def test_44_one_stored_set_is_used(self):
        make_subject(self.root, 's01', {'Artefact': {'NREM2': 1200.0, 'NREM3': 600.0}})
        dens, recs = self.ok()
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 30)
        self.assertEqual(recs['s01']['reject_types'], 'Artefact')

    def test_nothing_rejected_is_written_as_none(self):
        make_subject(self.root, 's01', {(0, 0): {'NREM2': 2400.0, 'NREM3': 1200.0}})
        dens, recs = self.ok()
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 60)
        self.assertEqual(recs['s01']['reject_types'], 'none')

    def test_several_sets_refused_without_the_option(self):
        make_subject(self.root, 's01', {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0},
                                        (1, 1): {'NREM2': 900.0, 'NREM3': 300.0}})
        self.refused(r'2 reject sets \(Arousal, Artefact; Artefact\).*--reject-types')

    def test_the_option_selects_a_stored_set(self):
        sets = {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0},
                (1, 1): {'NREM2': 900.0, 'NREM3': 300.0}}
        make_subject(self.root, 's01', sets)
        dens, recs = self.ok('--reject-types', ' Artefact , Arousal')   # order, spaces free
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 20)
        self.assertEqual(recs['s01']['reject_types'], 'Arousal,Artefact')
        dens, recs = self.ok('--reject-types', 'Artefact')
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 40)
        self.assertEqual(recs['s01']['reject_types'], 'Artefact')

    def test_a_named_set_not_stored_is_refused(self):
        make_subject(self.root, 's01', {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0}})
        self.refused(r'only for the reject set\(s\) Artefact, not for nothing rejected',
                     '--reject-types', '')
        self.refused(r'only for the reject set\(s\) Artefact, not for Arousal, Artefact, Move',
                     '--reject-types', 'Move,Artefact,Arousal')

    def test_none_names_the_set_that_rejects_nothing(self):
        # recordings.csv writes the empty set as 'none'; the option accepts it back.
        make_subject(self.root, 's01', {(0, 0): {'NREM2': 2400.0, 'NREM3': 1200.0},
                                        (1, 0): {'NREM2': 1800.0, 'NREM3': 600.0}})
        for name in ('none', 'None', ''):
            dens, recs = self.ok('--reject-types', name)
            self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 60)
            self.assertEqual(recs['s01']['reject_types'], 'none')

    def test_a_stored_key_not_in_turtlewave_spelling_is_refused(self):
        # One unsorted key: TurtleWave's density finds nothing under its own key, so the
        # bridge must not read it either.
        make_subject(self.root, 's01', {'Artefact,Arousal': {'NREM2': 900.0, 'NREM3': 300.0}})
        self.refused(r"'Artefact,Arousal'.*not in TurtleWave's key spelling")
        self.refused(r"not in TurtleWave's key spelling", '--reject-types', 'Arousal,Artefact')

    def test_two_spellings_of_one_set_are_refused_not_merged(self):
        make_subject(self.root, 's01', {'Arousal,Artefact': {'NREM2': 900.0, 'NREM3': 300.0},
                                        'Artefact,Arousal': {'NREM2': 60.0, 'NREM3': 60.0}})
        self.refused(r"'Artefact,Arousal'.*not in TurtleWave's key spelling")
        self.refused(r"not in TurtleWave's key spelling", '--reject-types', 'Arousal,Artefact')

    def test_an_empty_table_is_refused(self):
        make_subject(self.root, 's01', {(1, 0): {}})
        self.refused(r'stores no time under any reject set')

    def test_mixed_sets_across_databases_refused(self):
        make_subject(self.root, 's01', {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0}})
        make_subject(self.root, 's02', {'Arousal,Artefact': {'NREM2': 900.0, 'NREM3': 300.0}})
        self.refused(r'different reject sets.*\(Artefact: s01; Arousal, Artefact: s02\)')

    def test_a_set_named_and_stored_by_every_database_is_used_for_all(self):
        make_subject(self.root, 's01', {(1, 0): {'NREM2': 1800.0, 'NREM3': 600.0},
                                        (1, 1): {'NREM2': 900.0, 'NREM3': 300.0}})
        make_subject(self.root, 's02', {'Artefact': {'NREM2': 1200.0, 'NREM3': 600.0},
                                        'Arousal,Artefact,Move': {'NREM2': 60.0, 'NREM3': 60.0}})
        dens, recs = self.ok('--reject-types', 'Artefact')
        self.assertAlmostEqual(float(dens['s01']['E1']), 2 / 40)
        self.assertAlmostEqual(float(dens['s02']['E1']), 2 / 30)
        self.assertEqual({r['reject_types'] for r in recs.values()}, {'Artefact'})


if __name__ == '__main__':
    unittest.main()
