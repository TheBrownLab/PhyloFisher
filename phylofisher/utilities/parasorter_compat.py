#!/usr/bin/env python
import gzip
import os
import re
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


def to_tsv(args):
    '''
    Export ParaSorter2's saved (Parsed) decisions as <gene>_parsed.tsv files for apply_to_db.py.
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

    if 'concat_groups' in tables and conn.execute('SELECT 1 FROM concat_groups LIMIT 1').fetchone():
        print('Warning: this project has concatenation groups; PhyloFisher v1 cannot apply them and '
              'they are not exported.')

    os.makedirs(args.output, exist_ok=True)
    failed = []
    for gene in sorted(decisions):
        rows = decisions[gene]
        problems = []

        unrecognized = [tip for tip, _, _ in rows if forest_name(tip) is None]
        if unrecognized:
            problems.append(f'unrecognized tip names: {", ".join(unrecognized)}')

        # apply_to_db.py keeps only one ortholog per Unique ID, so more than one would silently drop data
        ortholog_ids = Counter(forest_name(tip).split('@')[-1].split('..')[0]
                               for tip, _, mode in rows if mode == 'o' and forest_name(tip))
        duplicated = sorted(uid for uid, count in ortholog_ids.items() if count > 1)
        if duplicated:
            problems.append(f'more than one ortholog for: {", ".join(duplicated)}')

        if problems:
            failed.append(gene)
            print(f'Error: {gene} not exported; ' + '; '.join(problems))
            continue

        untagged = tree_tips(tree_text(gene)) - {tip for tip, _, _ in rows}
        if untagged:
            print(f'Warning: {gene} has {len(untagged)} untagged tips, which apply_to_db.py will leave '
                  f'unchanged: {", ".join(sorted(untagged))}')

        with open(f'{args.output}/{gene}_parsed.tsv', 'w') as outfile:
            for tip, tax, mode in rows:
                outfile.write(f'{forest_name(tip)}\t{tax}\t{mode}\n')

    conn.close()

    not_reviewed = sorted(all_genes - set(decisions))
    print(f'Exported {len(decisions) - len(failed)} reviewed genes to {args.output}')
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
                   ready for apply_to_db.py -i <out_dir>.''')
    usage = ('parasorter_compat.py to_sqlite -i <sgt_constructor_out-local.tar.gz> [OPTIONS]\n'
             '       parasorter_compat.py to_tsv -i <project.sqlite> [OPTIONS]')

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
        in_help = 'Path to the ParaSorter2 project database (.sqlite)'
    else:
        parser.print_help()
        sys.exit('\nFirst argument must be to_sqlite or to_tsv')

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, out_dir=False, in_help=in_help)

    if mode == 'to_sqlite':
        to_sqlite(args)
    else:
        to_tsv(args)
