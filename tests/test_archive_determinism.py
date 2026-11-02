"""
Archives the aggregator writes are byte-identical when their contents are.

AmpliconRepository stores every version of a project, and can only share a
file between versions if it comes out byte-identical. Before this, every
archive recorded the time it was written (gzip header), each file's mtime,
uid/gid and owner names, and listed files in os.walk order, so an unchanged
CNVkit tarball or reconstruct log came out different on every run.

Run with: python -m pytest tests/
"""

import gzip
import os
import tarfile

from asa_aggregator import gzip_file_in_place, make_tarball

FILES = {
    'sample1.cns': b'chromosome\tstart\tend\n',
    'sub/sample1.cnr': b'cnr rows\n' * 50,
    'sub/deeper/sample1.call.cns': b'calls\n',
    'run.log': b'log line\n',
}


def _write_tree(root, order, mtime):
    for rel in order:
        path = os.path.join(root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'wb') as fh:
            fh.write(FILES[rel])
        os.utime(path, (mtime, mtime))
    return root


def _members(archive):
    with tarfile.open(archive, 'r:gz') as tar:
        return {m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}


def test_make_tarball_ignores_write_time_mtimes_and_walk_order(tmp_path, monkeypatch):
    a = _write_tree(str(tmp_path / 'a' / 'S1_cnvkit_output'), sorted(FILES), 1_000_000)
    b = _write_tree(str(tmp_path / 'b' / 'S1_cnvkit_output'), sorted(FILES, reverse=True), 2_000_000)
    os.chmod(os.path.join(b, 'run.log'), 0o600)  # a stricter umask elsewhere

    monkeypatch.setattr(gzip.time, 'time', lambda: 1_000_000)
    make_tarball(a, str(tmp_path / 'a.tar.gz'))
    monkeypatch.setattr(gzip.time, 'time', lambda: 2_000_000)
    make_tarball(b, str(tmp_path / 'b.tar.gz'))

    assert (tmp_path / 'a.tar.gz').read_bytes() == (tmp_path / 'b.tar.gz').read_bytes()


def test_make_tarball_keeps_contents_root_and_exclusions(tmp_path):
    src = _write_tree(str(tmp_path / 'cnvkit_output'), sorted(FILES), 1_000_000)
    make_tarball(src, str(tmp_path / 'out.tar.gz'), exclusions=('.log',),
                 root_name='S1_cnvkit_output')

    members = _members(str(tmp_path / 'out.tar.gz'))
    expected = {f'S1_cnvkit_output/{rel}': data for rel, data in FILES.items()
                if not rel.endswith('.log')}
    assert members == expected
    with tarfile.open(tmp_path / 'out.tar.gz', 'r:gz') as tar:
        assert [m.name for m in tar.getmembers()] == sorted(expected)


def test_gzip_file_in_place_ignores_write_time(tmp_path, monkeypatch):
    outputs = []
    for name, now in (('a', 1_000_000), ('b', 2_000_000)):
        monkeypatch.setattr(gzip.time, 'time', lambda now=now: now)
        d = tmp_path / name
        d.mkdir()
        path = d / 'S1_reconstruct.log'
        path.write_bytes(b'CoRAL log\n' * 100)
        assert gzip_file_in_place(str(path))
        assert not path.exists()
        outputs.append((d / 'S1_reconstruct.log.gz').read_bytes())

    assert outputs[0] == outputs[1]
    assert gzip.decompress(outputs[0]) == b'CoRAL log\n' * 100
    # The header still names the original file, as gzip.open did.
    assert b'S1_reconstruct.log\x00' in outputs[0][:64]


def test_archive_preserves_executability_and_unicode_paths(tmp_path):
    src = tmp_path / 'input'
    src.mkdir()
    for name, mode in [('analysis.sh', 0o700), ('résultats.tsv', 0o600)]:
        path = src / name
        path.write_bytes(b'contents\n')
        path.chmod(mode)
    archive = tmp_path / 'out.tar.gz'
    make_tarball(str(src), str(archive))
    with tarfile.open(archive) as tar:
        entries = {m.name: m for m in tar}
        assert entries['input/analysis.sh'].mode == 0o755
        assert entries['input/résultats.tsv'].mode == 0o644
        for member in entries.values():
            assert member.mtime == member.uid == member.gid == 0
            assert member.uname == member.gname == ''
            assert tar.extractfile(member).read() == b'contents\n'
    assert archive.read_bytes()[4:8] == b'\x00' * 4  # gzip mtime


def test_archive_omits_hidden_and_macos_files(tmp_path):
    src = tmp_path / 'input'
    for name in ['visible.txt', '.secret', '.hidden/data', '__MACOSX/data', 'sub/.secret']:
        path = src / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'contents\n')
    archive = tmp_path / 'out.tar.gz'
    make_tarball(str(src), str(archive))
    assert _members(archive) == {'input/visible.txt': b'contents\n'}
