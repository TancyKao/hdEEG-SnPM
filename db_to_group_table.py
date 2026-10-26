#!/usr/bin/env python3
"""
Bridge TurtleWave/wonambi per-subject output -> group-analysis tables.

For each subject folder under DATA (containing wonambi/neural_events.db and a
wonambi/*_clean_rebuilt.xml staging file), this, per event x frequency band:
  - selects the ONE detection run whose stored stage tokens lie inside the
    requested stages (--stages, default NREM2,NREM3). A joint run stored as
    'NREM2NREM3' and per-stage runs stored as 'NREM2' and 'NREM3' are both
    selected by NREM2,NREM3; a joint run is never split to answer a request for
    one of its stages. Two runs covering the same stage (a per-stage run and a
    joint run of the same night, or two methods) are REFUSED with an error,
    never pooled: their events largely duplicate each other;
  - takes per-channel COUNT, mean peak2peak amplitude and mean duration of the
    selected run's events;
  - computes density = count / analysed minutes, where analysed minutes is the
    time TurtleWave stored in its analysed_time table (TurtleWave >= 4.3) for the
    stages the selected run actually SEARCHED, intersected with --stages, under the
    reject set the detection stored. That set is read from the database, never
    assumed: the reject_types column (TurtleWave >= 4.4), else the reject_artifacts
    and reject_arousals flags of TurtleWave 4.3 ((1, 0) is 'Artefact', artefacts
    removed and arousals kept, which is what 4.3 detection stored; (1, 1) is
    'Arousal,Artefact'). A database storing one set is divided by that set's time. A
    database storing several is refused unless --reject-types names one, and a named
    set the database does not store is refused. Recordings whose densities would use
    different reject sets are refused, listing which recording used which set. A
    stored reject_types text that is not in TurtleWave's sorted key spelling was not
    written by TurtleWave and is refused. The searched stages are recovered as
    TurtleWave's own density code does: processing_status.stage for the run's method and band, else
    detection_runs.stages, joint tokens expanded to their stages. A run that
    searched only part of --stages (slow waves detected in NREM3 only, read with
    the default NREM2,NREM3) is divided by the time of the stages it searched, and
    a line "density for <band> over <stages> (run searched only these)" is printed
    for that recording. Only a database recording neither falls back, with a
    warning, to the stages its selected events are stored under (as TurtleWave
    and the reader do; --stages only when there is no event at all). Within one band every subject must be divided over the same stage
    set; a band whose subjects differ is refused, listing the subjects per set.
    A database without analysed_time (older TurtleWave) is REFUSED, as the
    hdEEG-SnPM Python reader and TurtleWave's event_density refuse it. Only with --allow-hypnogram-denominator does the
    bridge instead divide such a database by raw hypnogram minutes of the same
    stages, read from the first rater of the staging XML; that denominator is
    biased by each recording's artefact load, is announced with a warning per
    recording, and a run whose recordings would mix the two denominator types is
    refused. Which denominator was used, and over which stages, is printed per
    recording and band and written to recordings.csv; tables from separate runs
    must not be merged across denominator types or stage sets;
  - gives density 0 to a channel TurtleWave processed (processing_status) that
    produced no event, and sets every measure to missing on channels TurtleWave
    recorded as interpolated (detection_runs.interpolated_channels, TurtleWave
    >= 4.5), listing them in recordings.csv;
  - refuses, rather than writing a blank or zero row, a band whose events all
    lie under stage tokens --stages does not cover (a joint run cannot be split
    by stage), and a band with no event in scope whose processing_status lists
    more than one method (which method's channels are 0 is then undefined).
    A band with no event and no processed channel at all in a recording is
    written as a blank row, with a WARNING naming the subject and band (the
    Python reader stops there instead).
The run-selection and density rules follow the hdEEG-SnPM Python reader
(hdeeg_snpm.io.turtlewave) and TurtleWave's own event_density, re-implemented
here with the stdlib; the hypnogram fallback exists only in this bridge.

Output, in OUT:
  eventStat_<param>_<event_band>.csv   rows=subjects, cols=channels  (param in density/amplitude/duration)
  subjects.csv                         Subject, group, session
  channels.csv                         channel  (union across subjects)
  recordings.csv                       one row per subject x band: Subject, band,
                                       density_denominator (analysed_time|hypnogram),
                                       denominator_minutes, stages ('+'-joined),
                                       reject_types (the reject set of the analysed
                                       time, TurtleWave's sorted comma-joined key,
                                       'none' for the set that rejects nothing; blank
                                       for a hypnogram denominator),
                                       interpolated_channels (space-separated, per
                                       subject, repeated on each band row). The four
                                       denominator columns are blank for a band never
                                       detected in that subject.

Databases are opened read-only (SQLite URI mode=ro), so nothing is ever
written back into a TurtleWave database or its WAL file.

No dependency on the parameter/density CSVs (DB + XML only), so it survives the
planned TurtleWave outlier-removal / CSV-removal. Uses only the Python stdlib.

Usage:  python3 db_to_group_table.py <DATA_dir> <OUT_dir> [--group LABEL]
                                     [--stages NREM2,NREM3] [--bands LABEL,...]
                                     [--allow-hypnogram-denominator]
                                     [--reject-types SET]

--group is the cohort label written into subjects.csv. It defaults to 'group1',
which is only ever correct for a single-cohort run -- pass it explicitly when you
are building a table with more than one group, or every subject lands in the same
level and the group contrast has nothing to test.
--stages is the stage set analysed (comma-separated, or a joint token such as
NREM2NREM3). Default NREM2,NREM3. A run that searched only some of them is divided
by the time of the stages it searched (see above).
--bands restricts the output to some of the event-band labels below (default all).
--allow-hypnogram-denominator lets databases without analysed_time (TurtleWave
< 4.3) through, with density over raw hypnogram minutes. Off by default. Every
recording in the run must then lack analysed_time; a mix is refused.
--reject-types names the stored reject set whose analysed time divides density
("Artefact", "Artefact,Arousal"; order and spaces do not matter; "" or "none"
rejects nothing).
Needed only when a database stores more than one set; every database must store it.
"""
import os, sys, re, csv, sqlite3, glob, argparse
from pathlib import Path
import xml.etree.ElementTree as ET

EVENT_BANDS = [   # (event_type, freq_lower, freq_upper, label)
    ('spindle',   9.0, 12.0, 'spindle_9-12'),
    ('spindle',  12.0, 15.0, 'spindle_12-15'),
    ('slow_wave', 0.5, 1.25, 'sw_0.5-1.25'),
]
PARAMS = ['density', 'amplitude', 'duration']
DEFAULT_STAGES = ('NREM2', 'NREM3')
NREM_STAGES = set(DEFAULT_STAGES)              # kept for callers of stage_minutes()
STAGE_VOCAB = ('NREM1', 'NREM2', 'NREM3', 'Wake', 'REM')   # TurtleWave's stage vocabulary
HYPNO_FLAG = '--allow-hypnogram-denominator'


class BridgeError(Exception):
    """A database the bridge refuses to summarise; the message names the file."""


def warn(msg):
    print(f'WARNING: {msg}', file=sys.stderr)


def connect_ro(db_path):
    """Open a database strictly read-only (no WAL checkpoint, no journal written)."""
    return sqlite3.connect(Path(db_path).resolve().as_uri() + '?mode=ro', uri=True)


def stage_parts(stage):
    """['NREM2','NREM3'], 'NREM2' or the joint token 'NREM2NREM3' as a list of
    stages; a token outside the vocabulary ('Undefined') is one part."""
    if stage is None:
        return []
    if not isinstance(stage, str):
        return [str(s) for s in stage]
    parts, rest = [], stage
    while rest:
        hit = next((t for t in STAGE_VOCAB if rest.startswith(t)), None)
        if hit is None:
            return [stage]
        parts.append(hit)
        rest = rest[len(hit):]
    return parts


def parse_stages(text):
    """'NREM2,NREM3', 'NREM2 NREM3' or 'NREM2NREM3' -> ['NREM2', 'NREM3']."""
    out = []
    for item in re.split(r'[,\s+]+', text.strip()):
        for p in stage_parts(item) if item else []:
            if p not in STAGE_VOCAB:
                raise BridgeError(f"--stages: unknown stage '{p}' (known: {', '.join(STAGE_VOCAB)})")
            if p not in out:
                out.append(p)
    if not out:
        raise BridgeError('--stages is empty')
    return out


def stage_minutes(xml_path, stages=NREM_STAGES):
    """Raw hypnogram minutes in `stages`, from ONE rater (the first) of a Wonambi
    staging XML. Summing every rater would multiply the time by the rater count."""
    root = ET.parse(xml_path).getroot()
    raters = root.findall('.//rater')
    if len(raters) > 1:
        warn(f'{xml_path} holds {len(raters)} raters; the staging of the first '
             f"('{raters[0].get('name', '')}') is used.")
    scope = raters[0] if raters else root
    stages = set(stages)
    secs = 0.0
    for ep in scope.findall('.//epoch'):
        st = ep.find('stage')
        if st is not None and st.text in stages:
            secs += float(ep.find('epoch_end').text) - float(ep.find('epoch_start').text)
    return secs / 60.0


def _tables(con):
    return {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}


def _columns(con, table):
    return {r[1] for r in con.execute(f'PRAGMA table_info({table})')}


def interpolated_channels(con):
    """Channels TurtleWave recorded as interpolated (TurtleWave >= 4.5); else []."""
    if 'detection_runs' not in _tables(con) or \
            'interpolated_channels' not in _columns(con, 'detection_runs'):
        return []
    import json
    found = []
    for (cell,) in con.execute('SELECT interpolated_channels FROM detection_runs'):
        if cell:
            found += json.loads(cell)
    return list(dict.fromkeys(found))


def reject_key(reject_types):
    """A reject set as TurtleWave keys it (utils.reject_key): names stripped, sorted and
    joined by commas, 'Arousal,Artefact'; '' is the set that rejects nothing. Takes text
    ('Artefact, Arousal') or a list of names."""
    names = reject_types.split(',') if isinstance(reject_types, str) else list(reject_types)
    return ','.join(sorted({str(n).strip() for n in names if str(n).strip()}))


def reject_label(key):
    return key.replace(',', ', ') if key else 'nothing rejected'


def analysed_seconds(con, db_path, reject_types=None):
    """({stage: analysed seconds}, reject key) from TurtleWave's analysed_time, or None
    when the database predates it (no table).

    The reject set is the one the detection stored, never a library default: from the
    reject_types column (TurtleWave >= 4.4), else from the reject_artifacts /
    reject_arousals flags of TurtleWave 4.3 ((1, 0) is 'Artefact', (1, 1)
    'Arousal,Artefact', (0, 1) 'Arousal', (0, 0) '' = nothing rejected). One stored set
    is used as is; several are refused unless `reject_types` names one (--reject-types);
    a named set the database does not store is refused. A table with no usable set is
    refused too, and so is a stored reject_types text that is not TurtleWave's own key
    spelling (sorted, comma-joined, no spaces): TurtleWave did not write it, and
    TurtleWave's density would not find it either."""
    if 'analysed_time' not in _tables(con):
        return None
    subjects = [r[0] for r in con.execute('SELECT DISTINCT subject FROM analysed_time')]
    if len(subjects) > 1:
        raise BridgeError(f'{db_path}: analysed_time holds {len(subjects)} subjects; one '
                          'database per recording is expected.')
    cols = _columns(con, 'analysed_time')
    rows = []
    if subjects and 'reject_types' in cols:
        rows = [(k, st, sec) for k, st, sec in con.execute(
            'SELECT reject_types, stage, analysed_seconds FROM analysed_time WHERE subject = ?',
            subjects) if k is not None]
        odd = sorted({k for k, _, _ in rows if k != reject_key(k)})
        if odd:
            raise BridgeError(
                f'{db_path}: analysed_time stores reject set(s) '
                f"{'; '.join(repr(k) for k in odd)} that are not in TurtleWave's key "
                "spelling (names sorted, comma-joined, no spaces, e.g. 'Arousal,Artefact'), "
                'so they were written outside TurtleWave and TurtleWave\'s density would '
                'not find them. Fix the database; the bridge does not guess which set they '
                'mean.')
    elif subjects and {'reject_artifacts', 'reject_arousals'} <= cols:
        rows = [(reject_key([n for n, on in (('Artefact', a), ('Arousal', b)) if on]), st, sec)
                for a, b, st, sec in con.execute(
                    'SELECT reject_artifacts, reject_arousals, stage, analysed_seconds '
                    'FROM analysed_time WHERE subject = ?', subjects)]
    stored = sorted({k for k, _, _ in rows})
    listed = '; '.join(reject_label(k) for k in stored)
    if not stored:
        raise BridgeError(f'{db_path}: analysed_time stores no time under any reject set, so '
                          'there is nothing to divide the event count by and no density.')
    if reject_types is None:
        if len(stored) > 1:
            raise BridgeError(
                f'{db_path} stores analysed time for {len(stored)} reject sets ({listed}). '
                'Their densities cover different seconds; name the one to use with '
                '--reject-types (for example --reject-types "Artefact").')
        key = stored[0]
    else:
        key = '' if reject_types.strip().lower() == 'none' else reject_key(reject_types)
        if key not in stored:
            raise BridgeError(f'{db_path} stores analysed time only for the reject set(s) '
                              f'{listed}, not for {reject_label(key)} (--reject-types).')
    return {str(st): float(sec) for k, st, sec in rows if k == key}, key


def select_run(con, db_path, etype, flo, fhi, stages, label):
    """The stored stage tokens of the one run this band is read from.

    Returns (tokens, method): the covered stage tokens ([] only when the band has
    no event at all) and the run's method (None then). Raises BridgeError when two
    runs cover the request, or when the band's events all lie under tokens outside
    `stages`."""
    rows = con.execute(
        'SELECT method, stage, COUNT(*) FROM events WHERE event_type=? AND freq_lower=? '
        'AND freq_upper=? GROUP BY method, stage', (etype, flo, fhi)).fetchall()
    want = set(stages)
    scope = [r for r in rows if r[1] is not None and set(stage_parts(r[1])) <= want]
    other = [r for r in rows if r not in scope]
    methods = {m for m, _, _ in scope}
    parts = [p for _, t, _ in scope for p in stage_parts(t)]
    if len(methods) > 1 or len(parts) > len(set(parts)):
        runs = '; '.join(f"method {m}, stage token '{t}', {n} events" for m, t, n in scope)
        raise BridgeError(
            f'{db_path}: {label} ({etype} {flo:g}-{fhi:g} Hz) has more than one detection run '
            f"covering stages {'+'.join(stages)}: {runs}. They are not pooled, because the "
            'runs re-detect the same night and their events largely duplicate each other. '
            'Remove the unwanted run from the database, pick --stages so that only one run is '
            'covered, or leave this band out with --bands.')
    if other and scope:
        warn(f"{db_path}: {label}: {sum(n for _, _, n in other)} event(s) under stage token(s) "
             f"{sorted({str(t) for _, t, _ in other})} are outside --stages "
             f"{'+'.join(stages)} and are left out.")
    elif other and not scope:
        toks = sorted({str(t) for _, t, _ in other})
        wider = [t for t in toks if want & set(stage_parts(t))]
        why = (f"the stored token(s) {wider} also contain stages outside the request, and a "
               'joint run cannot be split by stage' if wider else
               'this band was not detected in the requested stages')
        raise BridgeError(
            f"{db_path}: {label} ({etype} {flo:g}-{fhi:g} Hz) has {sum(n for _, _, n in other)} "
            f"event(s), all under stage token(s) {toks}, none inside the requested stages "
            f"{'+'.join(stages)}: {why}. Writing this subject would give a blank row that "
            'looks like missing data. Pass --stages covering the stored token, or leave this '
            'band out with --bands.')
    return [t for _, t, _ in scope], (next(iter(methods)) if methods else None)


def processed_channels(con, db_path, etype, flo, fhi, method, label):
    """Channels TurtleWave processed successfully for this band (processing_status),
    so a channel that ran and found nothing gets density 0, not missing.

    `method` None means no event of the band selected a run; the roster is then
    usable only if processing_status lists a single method for the band."""
    need = {'channel', 'event_type', 'method', 'freq_lower', 'freq_upper', 'success'}
    if 'processing_status' not in _tables(con) or not need <= _columns(con, 'processing_status'):
        return set()
    rows = con.execute(
        'SELECT DISTINCT channel, method FROM processing_status WHERE success = 1 '
        'AND event_type = ? AND freq_lower = ? AND freq_upper = ?', (etype, flo, fhi)).fetchall()
    if method is None:
        listed = sorted({str(m) for _, m in rows})
        if len(listed) > 1:
            raise BridgeError(
                f'{db_path}: {label} ({etype} {flo:g}-{fhi:g} Hz) has no event in scope, and '
                f'processing_status lists {len(listed)} methods for it ({", ".join(listed)}), '
                'so which channels ran and found nothing (density 0) is undefined. Leave this '
                'band out with --bands, or read it with the hdEEG-SnPM Python reader naming '
                'method and stage.')
    return {str(ch) for ch, m in rows
            if method is None or str(m) == method or method in str(m).split('_')}


def _method_matches(stored, wanted):
    """TurtleWave's rule: processing_status/detection_runs store a run's method set
    '_'-joined; events.method holds one method. None matches anything."""
    return wanted is None or str(stored) == wanted or wanted in str(stored).split('_')


def searched_stages(con, etype, flo, fhi, method):
    """The stages the selected run actually SEARCHED, as TurtleWave's
    _run_stage_components recovers them: processing_status.stage for this method and
    band, else detection_runs.stages for this method and event type. Joint tokens are
    expanded to their component stages. A token that does not parse into known
    stages ('Undefined') is skipped, as TurtleWave skips it. [] when neither table
    can say."""
    import ast

    def strict(tok):
        parts = stage_parts(str(tok))
        return parts if all(p in STAGE_VOCAB for p in parts) else []

    found = []
    tabs = _tables(con)
    if 'processing_status' in tabs and {'stage', 'method', 'event_type', 'freq_lower',
                                        'freq_upper', 'success'} <= _columns(con, 'processing_status'):
        for m, tok in con.execute(
                'SELECT DISTINCT method, stage FROM processing_status WHERE success = 1 AND '
                'event_type = ? AND freq_lower = ? AND freq_upper = ?', (etype, flo, fhi)):
            if tok and _method_matches(m, method):
                found += strict(tok)
    if not found and 'detection_runs' in tabs and \
            {'stages', 'method', 'event_type'} <= _columns(con, 'detection_runs'):
        for m, rep in con.execute('SELECT method, stages FROM detection_runs WHERE event_type = ?',
                                  (etype,)):
            if not rep or not _method_matches(m, method):
                continue
            try:
                val = ast.literal_eval(rep)
            except (ValueError, SyntaxError):
                val = rep
            if isinstance(val, (list, tuple)):
                for v in val:
                    found += strict(v)
            elif isinstance(val, str):
                found += strict(val)
    return list(dict.fromkeys(found))


def subject_event_stats(db_path, xml_path, stages, bands, allow_hypnogram=False,
                        reject_types=None):
    """Return ({label: {channel: {density, amplitude, duration}}}, provenance dict).

    provenance['bands'][label] is (denominator type, minutes, [stages], note) for a
    band with a value in this recording, or None for a band never detected here;
    provenance['reject'] is the reject key of the analysed time divided by (None for a
    database without analysed_time)."""
    con = connect_ro(db_path)
    try:
        found = analysed_seconds(con, db_path, reject_types)
        denom, reject = found if found is not None else (None, None)
        if denom is None and not allow_hypnogram:
            raise BridgeError(
                f'{db_path} has no analysed_time table (it predates TurtleWave 4.3), so the '
                'artefact-free time density needs is not stored. Re-run detection with a '
                'current TurtleWave, or pass ' + HYPNO_FLAG + ' to divide by raw hypnogram '
                'minutes instead (biased by each recording\'s artefact load; every recording '
                'in the run must then lack analysed_time).')
        interp = interpolated_channels(con)
        out, prov = {}, {}
        for etype, flo, fhi, label in bands:
            tokens, method = select_run(con, db_path, etype, flo, fhi, stages, label)
            roster = processed_channels(con, db_path, etype, flo, fhi, method, label)
            rows = []
            if tokens:
                rows = con.execute(
                    'SELECT channel, COUNT(*), AVG(peak2peak_amp), AVG(duration) FROM events '
                    'WHERE event_type=? AND freq_lower=? AND freq_upper=? AND stage IN (%s) '
                    'GROUP BY channel' % ','.join('?' * len(tokens)),
                    (etype, flo, fhi, *tokens)).fetchall()
            if not rows and not roster:
                out[label], prov[label] = {}, None
                continue
            # Denominator stages: what the run SEARCHED, intersected with --stages. A
            # stage the run never searched cannot hold its events, so its time is not
            # part of the rate.
            token_parts = [p for t in tokens for p in stage_parts(t)]
            searched = searched_stages(con, etype, flo, fhi, method)
            note = ''
            if searched:
                dstages = [s for s in stages if s in searched]
                if not dstages:
                    raise BridgeError(
                        f"{db_path}: {label} ({etype} {flo:g}-{fhi:g} Hz): TurtleWave's "
                        f"processing record says the run searched only "
                        f"{'+'.join(searched)}, none of the requested stages "
                        f"{'+'.join(stages)}, so its density there is undefined. Pass "
                        '--stages covering the searched stages, or leave this band out '
                        'with --bands.')
                if len(dstages) < len(stages):
                    note = 'run searched only these'
            else:
                # No record of the searched stages (older database): the stages of the
                # stored event tokens stand in for them, as in TurtleWave and the reader.
                # --stages is used only when there is no token at all.
                dstages = [s for s in stages if s in token_parts] or list(stages)
                warn(f"{db_path}: {label}: the database records no searched stages for "
                     f"this run (no processing_status stage or detection_runs entry), so "
                     f"density is divided by the time of the stages its events are stored "
                     f"under: {'+'.join(dstages)}. A stage the run searched but found "
                     'nothing in leaves no event, so its time is missing from this '
                     'denominator and density may be over-estimated.')
            if denom is None:
                minutes = stage_minutes(xml_path, dstages)
                source = 'hypnogram'
            else:
                missing = [s for s in dstages if s not in denom]
                if missing:
                    raise BridgeError(
                        f"{db_path}: no analysed time stored under the reject set "
                        f"{reject_label(reject)} for stage(s) {', '.join(missing)} "
                        f"(stored: {', '.join(sorted(denom)) or 'none'}), "
                        f"so density for {label} over {'+'.join(dstages)} is undefined. Pass "
                        '--stages with the stages this recording was detected in. Raw '
                        'hypnogram time is not substituted on a database that stores '
                        'analysed time.')
                minutes = sum(denom[s] for s in dstages) / 60.0
                source = 'analysed_time'
            prov[label] = (source, round(minutes, 6), dstages, note)
            d = {ch: {'density': 0.0, 'amplitude': '', 'duration': ''} for ch in roster}
            for ch, n, amp, dur in rows:
                d[ch] = {'density': n / minutes if minutes > 0 else float('nan'),
                         'amplitude': amp, 'duration': dur}
            for ch in interp:
                if ch in d:
                    d[ch] = {p: '' for p in PARAMS}
            out[label] = d
    finally:
        con.close()
    return out, {'bands': prov, 'interpolated': interp, 'reject': reject}


def find_subjects(data_dir):
    subs = []
    for name in sorted(os.listdir(data_dir)):
        sdir = os.path.join(data_dir, name, 'wonambi')
        db = os.path.join(sdir, 'neural_events.db')
        xmls = sorted(glob.glob(os.path.join(sdir, '*_clean_rebuilt.xml')))
        if os.path.isfile(db) and xmls:
            subs.append((name, db, xmls[0]))
    return subs


def natkey(ch):  # E1,E2,...,E10,...,Cz  natural order
    m = re.match(r'E(\d+)$', ch)
    return (0, int(m.group(1))) if m else (1, ch)


def run(a):
    stages = parse_stages(a.stages)
    bands = EVENT_BANDS
    if a.bands:
        pick = [b for b in re.split(r'[,\s]+', a.bands.strip()) if b]
        unknown = sorted(set(pick) - {lab for *_, lab in EVENT_BANDS})
        if unknown:
            raise BridgeError(f"--bands: unknown label(s) {unknown}; known: "
                              f"{[lab for *_, lab in EVENT_BANDS]}")
        bands = [b for b in EVENT_BANDS if b[3] in pick]
    os.makedirs(a.out_dir, exist_ok=True)

    subs = find_subjects(a.data_dir)
    if not subs:
        raise BridgeError(f'No subjects (wonambi/neural_events.db + xml) under {a.data_dir}')
    print(f'Subjects: {[s[0] for s in subs]}')
    print(f"Stages: {'+'.join(stages)}")

    # stats[label][param][subject][channel] = value ; collect channel union
    stats = {lab: {p: {} for p in PARAMS} for *_, lab in bands}
    channels = set()
    meta, recs = [], []
    for name, db, xml in subs:
        ev, prov = subject_event_stats(db, xml, stages, bands, a.allow_hypnogram_denominator,
                                       a.reject_types)
        if any(v and v[0] == 'hypnogram' for v in prov['bands'].values()):
            warn(f'{name}: database predates TurtleWave 4.3 (no analysed_time), so density '
                 'is divided by raw hypnogram minutes, which include artefact time; density '
                 'is under-estimated in proportion to this recording\'s artefact load.')
        parts = []
        for *_, lab in bands:
            v = prov['bands'][lab]
            if v is None:
                warn(f'{name}: {lab} was never detected in this recording (no event and no '
                     'processed-channel record for it), so its row is blank in every '
                     f'eventStat_*_{lab}.csv table.')
                parts.append(f'{lab}: not detected')
                continue
            kind = (f"TurtleWave analysed time ({reject_label(prov['reject'])})"
                    if v[0] == 'analysed_time' else 'hypnogram')
            parts.append(f"{lab}: {len(ev[lab])} channels, density over {kind} "
                         f"{'+'.join(v[2])} = {v[1]:.1f} min")
            if v[3]:
                print(f"  {name}: density for {lab} over {'+'.join(v[2])} ({v[3]})")
        print(f'  {name}: ' + '; '.join(parts)
              + (f"; interpolated (set missing): {' '.join(prov['interpolated'])}"
                 if prov['interpolated'] else ''))
        for lab in ev:
            for ch, vals in ev[lab].items():
                channels.add(ch)
                for p in PARAMS:
                    stats[lab][p].setdefault(name, {})[ch] = vals[p]
        sess = name.split('_')[-1] if '_' in name else ''
        meta.append((name, a.group, sess))
        interp = ' '.join(prov['interpolated'])
        for *_, lab in bands:
            v = prov['bands'][lab]
            used = v is not None and v[0] == 'analysed_time'
            recs.append((name, lab) + ((v[0], v[1], '+'.join(v[2])) if v else ('', '', ''))
                        + ((prov['reject'] or 'none') if used else '', interp))

    kinds = {}
    for name, _, kind, *_ in recs:
        if kind and name not in kinds.get(kind, []):
            kinds.setdefault(kind, []).append(name)
    if len(kinds) > 1:
        raise BridgeError(
            'recordings resolve to different density denominators, which would put two '
            'scales in one table: ' + '; '.join(f'{k}: {", ".join(v)}' for k, v in
                                                 sorted(kinds.items()))
            + '. Run the two kinds separately and do not compare them, or re-run detection '
            'on the old databases with a current TurtleWave. Nothing was written.')
    rsets = {}
    for name, _, kind, _, _, key, _ in recs:
        if kind == 'analysed_time' and name not in rsets.get(key, []):
            rsets.setdefault(key, []).append(name)
    if len(rsets) > 1:
        raise BridgeError(
            'The databases were detected with different reject sets, so their densities are '
            'divided by time cleaned in different ways ('
            + '; '.join(f'{reject_label(k)}: {", ".join(v)}' for k, v in rsets.items())
            + ') and are not comparable. Name one set every database stores with '
            '--reject-types, or redetect so the runs match. Nothing was written.')
    for *_, lab in bands:
        sets = {}
        for name, b, kind, _, st, _, _ in recs:
            if b == lab and kind:
                sets.setdefault(st, []).append(name)
        if len(sets) > 1:
            raise BridgeError(
                f'{lab}: subjects\' densities are taken over different stage sets, which '
                'would put rates over different parts of the night in one table: '
                + '; '.join(f'{st}: {", ".join(v)}' for st, v in sorted(sets.items()))
                + '. Pass --stages so that every subject is read over the same stages, '
                're-run detection so the runs search the same stages, or leave this band '
                'out with --bands. Nothing was written.')

    chans = sorted(channels, key=natkey)
    subjects = [s[0] for s in subs]

    # write per (param, event-band) subject x channel matrices
    for *_, lab in bands:
        for p in PARAMS:
            path = os.path.join(a.out_dir, f'eventStat_{p}_{lab}.csv')
            with open(path, 'w', newline='') as f:
                w = csv.writer(f); w.writerow(['Subject'] + chans)
                for s in subjects:
                    row = [s] + [stats[lab][p].get(s, {}).get(ch, '') for ch in chans]
                    w.writerow(row)
    # subjects.csv + channels.csv + recordings.csv
    with open(os.path.join(a.out_dir, 'subjects.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['Subject', 'group', 'session']); w.writerows(meta)
    with open(os.path.join(a.out_dir, 'channels.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['channel']); w.writerows([[c] for c in chans])
    with open(os.path.join(a.out_dir, 'recordings.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Subject', 'band', 'density_denominator', 'denominator_minutes',
                    'stages', 'reject_types', 'interpolated_channels'])
        w.writerows(recs)

    print(f'Wrote {len(bands)*len(PARAMS)} matrices + subjects.csv + channels.csv + '
          f'recordings.csv ({len(subjects)} subjects x {len(chans)} channels) to {a.out_dir}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data_dir'); ap.add_argument('out_dir')
    ap.add_argument('--group', default='group1',
                    help="cohort label written into subjects.csv (default: group1)")
    ap.add_argument('--stages', default=','.join(DEFAULT_STAGES),
                    help='stages analysed, e.g. NREM2,NREM3 (default) or NREM3')
    ap.add_argument('--bands', default='',
                    help='comma-separated event-band labels to write (default: all)')
    ap.add_argument(HYPNO_FLAG, action='store_true',
                    help='accept databases without analysed_time (TurtleWave < 4.3) and divide '
                         'their density by raw hypnogram minutes (biased; off by default)')
    ap.add_argument('--reject-types', default=None, metavar='SET',
                    help='the stored reject set whose analysed time divides density, e.g. '
                         '"Artefact" or "Artefact,Arousal" ("" or "none" rejects nothing, '
                         'as recordings.csv writes it). Needed only '
                         'when a database stores several sets (default: the one set stored)')
    a = ap.parse_args()
    try:
        run(a)
    except BridgeError as e:
        sys.exit(f'ERROR: {e}')


if __name__ == '__main__':
    main()
