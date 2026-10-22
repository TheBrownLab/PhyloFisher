#!/usr/bin/env python
import gzip
import os
import re
import shutil
import sqlite3
import sys
import tarfile
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

from phylofisher import help_formatter

TREE_EXTS = ('raxml.support', 'treefile')
ALN_EXTS = ('final', 'trimmed', 'aln')
GENE_FILE_RE = re.compile(r'^(.+)\.(raxml\.support|treefile|final|trimmed|aln)$')
SIDE_FILES = ('metadata.tsv', 'input_metadata.tsv', 'tree_colors.tsv')
LEAF_RE = re.compile(r'[(,]\s*([^(),:;\s]+)')


def collect_from_tarball(tar_path):
    '''
    Collect project files from a v1 sgt_constructor_out_<date>-local.tar.gz.

    :param tar_path: path to the tarball
    :type tar_path: str
    :return: dictionary of project-relative path -> file bytes
    :rtype: dict
    '''
    files = {}
    with tarfile.open(tar_path, 'r:gz') as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            parts = Path(member.name).parts
            # Members live under <out_dir>-local/, so drop that leading directory
            if len(parts) == 2 and parts[1] in SIDE_FILES:
                rel_path = parts[1]
            elif len(parts) == 3 and parts[1] == 'trees' and GENE_FILE_RE.match(parts[2]):
                rel_path = f'trees/{parts[2]}'
            else:
                continue
            files[rel_path] = tar.extractfile(member).read()
    return files


def default_sqlite_path(in_path):
    '''
    Mirror PhyloFisher v2's naming: <out_dir>-local.sqlite next to the sgt_constructor output.
    '''
    base = os.path.abspath(in_path).rstrip('/')
    for suffix in ('.tar.gz', '.tgz'):
        if base.endswith(suffix):
            base = base[:-len(suffix)]
    if not base.endswith('-local'):
        base = f'{base}-local'
    return f'{base}.sqlite'


def to_sqlite(args):
    '''
    Write a ParaSorter2 project database from v1 sgt_constructor output.
    '''
    if not (os.path.isfile(args.input) and tarfile.is_tarfile(args.input)):
        sys.exit(f'{args.input} is not an sgt_constructor -local.tar.gz')
    files = collect_from_tarball(args.input)

    genes = {GENE_FILE_RE.match(p.split('/')[-1]).group(1) for p in files if p.startswith('trees/')}
    tree_genes = {g for g in genes if any(f'trees/{g}.{ext}' in files for ext in TREE_EXTS)}
    if not tree_genes:
        sys.exit(f'No gene trees (*.raxml.support or *.treefile) found in {args.input}')
    for gene in sorted(genes - tree_genes):
        print(f'Warning: {gene} has no tree and will not be listed in ParaSorter')
    for name in SIDE_FILES:
        if name not in files:
            print(f'Warning: {name} not found in the tarball; ParaSorter will show raw tip names for these taxa.')

    out_path = args.output or default_sqlite_path(args.input)
    if os.path.exists(out_path):
        if not args.force:
            sys.exit(f'{out_path} already exists. Use --force to overwrite it '
                     '(this discards any ParaSorter decisions saved in it).')
        os.remove(out_path)

    conn = sqlite3.connect(out_path)
    with conn:
        conn.execute('CREATE TABLE files (path TEXT PRIMARY KEY, content BLOB NOT NULL)')
        conn.executemany('INSERT INTO files (path, content) VALUES (?, ?)',
                         ((rel, gzip.compress(content)) for rel, content in sorted(files.items())))
    conn.close()

    print(f'Wrote {len(tree_genes)} gene trees ({len(files)} files) to {out_path}')


def tree_tips(newick):
    '''
    Return the leaf names of a Newick tree.
    '''
    return set(LEAF_RE.findall(newick))


def forest_name(tip):
    '''
    Convert a raw tip name into the forest.py-style name apply_to_db.py parses.

    apply_to_db.py takes the Unique ID from after the last "@", treats names with exactly
    three underscores as new fisher candidates with their quality (e.g. q1c) in field 2,
    and treats names containing ".." as database paralogs.

    :param tip: raw tip name from the gene tree
    :type tip: str
    :return: name for the _parsed.tsv, or None if the tip is not in a recognized format
    :rtype: str
    '''
    n_underscores = tip.count('_')
    if n_underscores == 3:
        # New candidate: <UniqueID>_<n>_<SBH|BBH|HMM>_q<N><c|n>
        unique_id, _, route, quality = tip.split('_')
        return f'{unique_id}_{route}_{quality}_NA@{unique_id}'
    if n_underscores == 0:
        # Database ortholog (<UniqueID>) or paralog (<UniqueID>..p<N>)
        return tip
    return None


def concat_name(unique_id, group_id):
    '''
    Names for a concatenation group's merged sequence: (record name in the fisher <gene>.fas,
    name for the _parsed.tsv).

    Both have three underscores, so apply_to_db.py treats the merged sequence as a new candidate
    for unique_id; it matches the two on <UniqueID>_<field 2 of the tsv name>, so "merged<N>" has
    to be unique within the gene and must not look like a fisher quality (q1c, ...).
    '''
    return (f'{unique_id}_{group_id}_CONCAT_merged{group_id}',
            f'{unique_id}_CONCAT_merged{group_id}_NA@{unique_id}')


def read_concat_groups(conn, tables):
    '''
    Read ParaSorter2's concatenation groups and their merged sequences.

    :return: dictionary of gene -> list of (group_id, member tips in priority order, mode,
             merged sequence or None if ParaSorter has not saved an up-to-date one)
    :rtype: dict
    '''
    if 'concat_groups' not in tables:
        return {}
    members = defaultdict(lambda: defaultdict(list))
    modes = {}
    # rowid order is the order ParaSorter wrote each group's tips in, i.e. merge priority
    for gene, group_id, tip, mode in conn.execute(
            'SELECT gene_id, group_id, tip, mode FROM concat_groups ORDER BY gene_id, group_id, rowid'):
        members[gene][group_id].append(tip)
        modes[(gene, group_id)] = mode

    merged = {}
    if 'concat_sequences' in tables:
        for gene, group_id, tips, sequence in conn.execute(
                'SELECT gene_id, group_id, tips, sequence FROM concat_sequences'):
            merged[(gene, group_id)] = (tips.split('\n'), sequence)

    groups = defaultdict(list)
    for gene, gene_groups in members.items():
        for group_id, tips in sorted(gene_groups.items()):
            saved_tips, sequence = merged.get((gene, group_id), (None, None))
            groups[gene].append((group_id, tips, modes[(gene, group_id)],
                                 sequence if saved_tips == tips else None))
    return groups


def copy_fisher_dir(fisher_dir, fisher_out):
    '''
    Copy the top-level files of a fisher.py output directory (the <gene>.fas files apply_to_db.py
    reads, plus its small tables); tmp/ is not needed by apply_to_db.py and is skipped.
    '''
    if not os.path.isdir(fisher_dir):
        sys.exit(f'{fisher_dir} is not a fisher.py output directory')
    if os.path.abspath(fisher_dir) == os.path.abspath(fisher_out):
        sys.exit('--fisher_out must be different from --fisher_dir')
    os.makedirs(fisher_out, exist_ok=True)
    for entry in os.scandir(fisher_dir):
        if entry.is_file():
            shutil.copy2(entry.path, fisher_out)


def to_tsv(args):
    '''
    Export ParaSorter2's saved (Parsed) decisions as <gene>_parsed.tsv files for apply_to_db.py.

    Each concatenation group becomes one new candidate: its merged sequence is appended to the
    gene's .fas in a copy of the fisher.py output (--fisher_out), tagged with the group's o/p, and
    its member tips are tagged d.
    '''
    if not os.path.isfile(args.input):
        sys.exit(f'{args.input} does not exist')
    conn = sqlite3.connect(f'file:{os.path.abspath(args.input)}?mode=ro', uri=True)
    tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
    if 'files' not in tables:
        sys.exit(f'{args.input} is not a ParaSorter project database (no "files" table)')
    if 'tagging' not in tables:
        sys.exit(f'{args.input} has no saved ParaSorter decisions (no "tagging" table)')

    decisions = defaultdict(list)
    for gene, tip, tax, mode in conn.execute(
            "SELECT gene_id, tip, class, mode FROM tagging WHERE kind = 'parsed' ORDER BY gene_id, tip"):
        decisions[gene].append((tip, tax, mode))

    all_genes = set()
    for (path,) in conn.execute("SELECT path FROM files WHERE path LIKE 'trees/%'"):
        m = GENE_FILE_RE.match(path[len('trees/'):])
        if m and m.group(2) in ('raxml.support', 'treefile'):
            all_genes.add(m.group(1))

    def tree_text(gene):
        for ext in TREE_EXTS:
            row = conn.execute('SELECT content FROM files WHERE path = ?', (f'trees/{gene}.{ext}',)).fetchone()
            if row:
                return gzip.decompress(row[0]).decode()
        return ''

    concat_groups = read_concat_groups(conn, tables)
    if concat_groups and not args.fisher_dir:
        sys.exit('This project has concatenation groups. Re-run with -fi <fisher_dir> so their merged '
                 'sequences can be added to a copy of the fisher.py output.')
    fisher_out = args.fisher_out or f'{args.output}/fisher_out'
    if args.fisher_dir:
        copy_fisher_dir(args.fisher_dir, fisher_out)

    os.makedirs(args.output, exist_ok=True)
    genes = sorted(set(decisions) | set(concat_groups))
    failed = []
    for gene in genes:
        groups = concat_groups.get(gene, [])
        grouped_tips = {tip for _, tips, _, _ in groups for tip in tips}
        # ParaSorter pauses a tip's own o/p/d while it is grouped, but a Parsed save from before it
        # was grouped can still be there; the group decides what happens to it
        rows = [row for row in decisions[gene] if row[0] not in grouped_tips]
        problems = []

        unrecognized = [tip for tip, _, _ in rows if forest_name(tip) is None]
        unrecognized += [tip for tip in sorted(grouped_tips) if forest_name(tip) is None]
        if unrecognized:
            problems.append(f'unrecognized tip names: {", ".join(unrecognized)}')
            failed.append(gene)
            print(f'Error: {gene} not exported; ' + '; '.join(problems))
            continue

        out_rows = [(forest_name(tip), tax, mode) for tip, tax, mode in rows]
        merged_records = []
        for group_id, tips, mode, sequence in groups:
            unique_ids = {forest_name(tip).split('@')[-1].split('..')[0] for tip in tips}
            if len(unique_ids) > 1:
                problems.append(f'concatenation group {group_id} mixes taxa: {", ".join(sorted(unique_ids))}')
                continue
            if not sequence:
                problems.append(f'no up-to-date merged sequence for concatenation group {group_id} '
                                f'(open {gene} in ParaSorter to save it; the gene needs an untrimmed .aln '
                                f'alignment containing every member)')
                continue
            unique_id = unique_ids.pop()
            for tip in tips:
                if '..' in tip:
                    print(f'Warning: {gene}: {tip} is a database paralog in concatenation group {group_id}; '
                          f'apply_to_db.py cannot remove paralogs, so it will stay in the database too')
                out_rows.append((forest_name(tip), '', 'd'))
            record_name, tsv_name = concat_name(unique_id, group_id)
            out_rows.append((tsv_name, '', mode))
            merged_records.append((record_name, sequence))

        # apply_to_db.py keeps only one ortholog per Unique ID, so more than one would silently drop data
        ortholog_ids = Counter(name.split('@')[-1].split('..')[0] for name, _, mode in out_rows if mode == 'o')
        duplicated = sorted(uid for uid, count in ortholog_ids.items() if count > 1)
        if duplicated:
            problems.append(f'more than one ortholog for: {", ".join(duplicated)}')

        fas_path = f'{fisher_out}/{gene}.fas'
        if merged_records and not os.path.isfile(fas_path):
            problems.append(f'{gene}.fas not found in {args.fisher_dir}')

        if problems:
            failed.append(gene)
            print(f'Error: {gene} not exported; ' + '; '.join(problems))
            continue

        untagged = tree_tips(tree_text(gene)) - {tip for tip, _, _ in rows} - grouped_tips
        if untagged:
            print(f'Warning: {gene} has {len(untagged)} untagged tips, which apply_to_db.py will leave '
                  f'unchanged: {", ".join(sorted(untagged))}')

        if merged_records:
            with open(fas_path, 'a') as outfile:
                for record_name, sequence in merged_records:
                    outfile.write(f'>{record_name}\n{sequence}\n')

        with open(f'{args.output}/{gene}_parsed.tsv', 'w') as outfile:
            for name, tax, mode in out_rows:
                outfile.write(f'{name}\t{tax}\t{mode}\n')

    conn.close()

    not_reviewed = sorted(all_genes - set(genes))
    print(f'Exported {len(genes) - len(failed)} reviewed genes to {args.output}')
    if args.fisher_dir:
        print(f'Wrote fisher.py output with merged sequences to {fisher_out}; '
              f'run apply_to_db.py -i {args.output} -fi {fisher_out}')
    if not_reviewed:
        print(f'{len(not_reviewed)} genes have no saved (Parsed) decisions and were skipped')
    if failed:
        sys.exit(f'{len(failed)} genes were not exported. Fix them in ParaSorter and re-run.')


if __name__ == '__main__':
    description = textwrap.dedent('''\
        Backwards compatibility between PhyloFisher v1 and ParaSorter2.

        to_sqlite: package sgt_constructor.py output into a ParaSorter2 project database
                   (<sgt_constructor_out>-local.sqlite).
        to_tsv:    export the decisions saved in ParaSorter2 as <gene>_parsed.tsv files,
                   ready for apply_to_db.py -i <out_dir>. With -fi, also writes a copy of the
                   fisher.py output with concatenation groups' merged sequences added, for
                   apply_to_db.py -fi.''')
    usage = ('parasorter_compat.py to_sqlite -i <sgt_constructor_out-local.tar.gz> [OPTIONS]\n'
             '       parasorter_compat.py to_tsv -i <project.sqlite> [-fi <fisher_dir>] [OPTIONS]')

    mode = sys.argv.pop(1) if len(sys.argv) > 1 and sys.argv[1] in ('to_sqlite', 'to_tsv') else None
    name = f'parasorter_compat.py {mode}' if mode else 'parasorter_compat.py'
    parser, optional, required = help_formatter.initialize_argparse(name=name, desc=description, usage=usage)

    if mode == 'to_sqlite':
        optional.add_argument('-o', '--output', type=str, metavar='<out.sqlite>', default=None,
                              help=textwrap.dedent("""\
                              Path to the project database to write.
                              Default: <input>-local.sqlite next to the input
                              """))
        optional.add_argument('--force', action='store_true',
                              help='Overwrite an existing project database')
        in_help = 'Path to the sgt_constructor.py output -local.tar.gz'
    elif mode == 'to_tsv':
        today_date = help_formatter.today.strftime("%b.%d.%Y")
        optional.add_argument('-o', '--output', type=str, metavar='<out_dir>', default=f'parasorter_out_{today_date}',
                              help=textwrap.dedent("""\
                              Directory to write <gene>_parsed.tsv files to.
                              Default: ./parasorter_out_<M.D.Y>
                              """))
        optional.add_argument('-fi', '--fisher_dir', type=str, metavar='<fisher_dir>', default=None,
                              help=textwrap.dedent("""\
                              Path to the fisher.py output directory the project was built from.
                              Required if the project has concatenation groups.
                              """))
        optional.add_argument('-fo', '--fisher_out', type=str, metavar='<fisher_out>', default=None,
                              help=textwrap.dedent("""\
                              Where to write the copy of <fisher_dir> with merged sequences added,
                              to give to apply_to_db.py -fi.
                              Default: <out_dir>/fisher_out
                              """))
        in_help = 'Path to the ParaSorter2 project database (.sqlite)'
    else:
        parser.print_help()
        sys.exit('\nFirst argument must be to_sqlite or to_tsv')

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, out_dir=False, in_help=in_help)

    if mode == 'to_sqlite':
        to_sqlite(args)
    else:
        to_tsv(args)
