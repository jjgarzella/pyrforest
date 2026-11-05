"""Build and exercise a separate consumer using a built pyrforest wheel.

Usage: sage -python tests/build_native_pn_consumer.py path/to/pyrforest.whl
"""

import hashlib
import importlib
import importlib.machinery
import os
import subprocess
import sys
import tempfile
import textwrap
import zipfile
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Distribution, Extension
from setuptools.command.build_ext import build_ext


def _sage_includes():
    try:
        from sage.config import get_include_dirs

        return [str(path) for path in get_include_dirs()]
    except ImportError:
        from sage.env import sage_include_directories

        return [str(path) for path in sage_include_directories()]


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_inside(path, root, description):
    if path is None:
        raise AssertionError(f"{description} was not loaded")
    path = Path(path).resolve()
    root = Path(root).resolve()
    if not path.is_file():
        raise AssertionError(f"{description} is not a file: {path}")
    try:
        path.relative_to(root)
    except ValueError as error:
        raise AssertionError(
            f"{description} resolved outside the expected wheel directory: {path}"
        ) from error
    return path


def _validated_producer_paths(package_path, extension_path, package_dir):
    package_path = _file_inside(package_path, package_dir, "producer package")
    extension_path = _file_inside(
        extension_path, package_dir, "compiled producer extension"
    )
    if not any(
        extension_path.name.endswith(suffix)
        for suffix in importlib.machinery.EXTENSION_SUFFIXES
    ):
        raise AssertionError(
            f"producer did not load from a compiled extension: {extension_path}"
        )
    return package_path, extension_path


def build_consumer(wheel_path):
    wheel_path = Path(wheel_path).resolve()
    source = Path(__file__).with_name("native_pn_consumer.pyx").resolve()
    if not wheel_path.is_file():
        raise FileNotFoundError(wheel_path)
    print(f"WHEEL_ARCHIVE {wheel_path} sha256={_sha256(wheel_path)}")
    with tempfile.TemporaryDirectory(prefix="pyrforest-pn-consumer-") as work:
        root = Path(work)
        unpacked = root / "wheel"
        unpacked.mkdir()
        with zipfile.ZipFile(wheel_path) as archive:
            archive.extractall(unpacked)
        package = unpacked / "pyrforest"
        header = package / "native_pn.h"
        pxd = package / "native_pn.pxd"
        if not header.is_file() or not pxd.is_file():
            raise AssertionError("wheel omitted the PN GMP header or Cython declaration")

        extension = Extension(
            "native_pn_consumer",
            [str(source)],
            include_dirs=[str(package), *_sage_includes()],
            libraries=["gmp"],
            extra_compile_args=["-std=gnu11"],
        )
        extensions = cythonize(
            [extension],
            include_path=[str(unpacked), *_sage_includes()],
            build_dir=str(root / "cython"),
            compiler_directives={"language_level": 3},
            quiet=True,
        )
        distribution = Distribution({"name": "pyrforest-pn-consumer-test",
                                     "ext_modules": extensions})
        command = build_ext(distribution)
        command.build_lib = str(root / "lib")
        command.build_temp = str(root / "build")
        command.ensure_finalized()
        command.run()

        # A conflicting package on PYTHONPATH reproduces the original failure
        # mode. The isolated child must ignore it and prove that the extracted
        # wheel supplied both the package and compiled producer extension.
        ambient = root / "ambient"
        ambient_package = ambient / "pyrforest"
        ambient_package.mkdir(parents=True)
        marker = root / "ambient-producer-imported"
        (ambient_package / "__init__.py").write_text(
            textwrap.dedent(f"""\
            from pathlib import Path
            Path({str(marker)!r}).write_text("imported")
            raise RuntimeError("ambient pyrforest producer sentinel imported")
            """)
        )
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            value for value in (env.get("PYTHONPATH"), str(ambient)) if value
        )
        runner = subprocess.run(
            [
                sys.executable,
                "-I",
                str(Path(__file__).resolve()),
                "--isolated-consumer",
                str(unpacked),
                str(command.build_lib),
            ],
            cwd=root,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if runner.stdout:
            print(runner.stdout, end="")
        if runner.returncode:
            if runner.stderr:
                print(runner.stderr, file=sys.stderr, end="")
            raise RuntimeError(
                f"isolated wheel consumer exited with status {runner.returncode}"
            )
        if marker.exists():
            raise AssertionError("ambient producer sentinel was imported")
        print(f"AMBIENT_SENTINEL {ambient_package} ignored under PYTHONPATH")


def _run_isolated_consumer(wheel_root, consumer_lib):
    wheel_root = Path(wheel_root).resolve()
    consumer_lib = Path(consumer_lib).resolve()
    package_dir = wheel_root / "pyrforest"
    header = package_dir / "native_pn.h"
    pxd = package_dir / "native_pn.pxd"
    for required in (header, pxd):
        if not required.is_file():
            raise AssertionError(f"wheel omitted required ABI file: {required}")

    # The extension uses Sage's Cython and symbolic runtime; initialize Sage
    # before importing the producer package, as documented for consumers.
    from sage.all import Integers, Matrix, PolynomialRing, ZZ

    # The child runs from a temporary working directory under Python's isolated
    # mode. Only these explicit paths can supply the producer and consumer.
    sys.path.insert(0, str(wheel_root))
    sys.path.insert(0, str(consumer_lib))
    package = importlib.import_module("pyrforest")
    extension = importlib.import_module("pyrforest.rforest")
    consumer = importlib.import_module("native_pn_consumer")
    package_path, extension_path = _validated_producer_paths(
        package.__file__, extension.__file__, package_dir
    )
    package_dir_loaded = Path(next(iter(package.__path__))).resolve()
    if package_dir_loaded != package_dir.resolve():
        raise AssertionError(
            f"producer package directory came from {package_dir_loaded}, "
            f"expected extracted wheel at {package_dir.resolve()}"
        )
    consumer_path = _file_inside(
        consumer.__file__, consumer_lib, "compiled consumer extension"
    )
    if not any(
        consumer_path.name.endswith(suffix)
        for suffix in importlib.machinery.EXTENSION_SUFFIXES
    ):
        raise AssertionError(
            f"consumer was not loaded from a compiled extension: {consumer_path}"
        )

    print(f"PRODUCER_PACKAGE {package_path} sha256={_sha256(package_path)}")
    print(f"PRODUCER_EXTENSION {extension_path} sha256={_sha256(extension_path)}")
    print(f"WHEEL_HEADER {header.resolve()} sha256={_sha256(header)}")
    print(f"WHEEL_PXD {pxd.resolve()} sha256={_sha256(pxd)}")
    consumer_source = Path(__file__).with_name("native_pn_consumer.pyx").resolve()
    print(f"CONSUMER_SOURCE {consumer_source} sha256={_sha256(consumer_source)}")
    print(f"CONSUMER_EXTENSION {consumer_path} sha256={_sha256(consumer_path)}")

    from pyrforest import remainder_forest_pn, remainder_forest_pn_native

    ring = PolynomialRing(ZZ, names=("P", "x"))
    P, x = ring.gens()
    matrix = Matrix(ring, [[1 + P + x, -2 + 3 * P], [4 - P * x, 2 + x]])
    initial = Matrix(
        ring,
        [[2**360 + 2 + (2**330 + 1) * P,
          -(2**315 + 3) + 2 * P]],
    )
    moduli = [2**501 + 17, 2**431 + 93]
    endpoints = [4, 6]
    indices = ("first", "second")
    noncommuting_ring = PolynomialRing(Integers(101), "P")
    transition_one = Matrix(
        noncommuting_ring,
        matrix.nrows(),
        matrix.ncols(),
        [noncommuting_ring(entry.subs({x: 1})) for entry in matrix.list()],
    )
    transition_two = Matrix(
        noncommuting_ring,
        matrix.nrows(),
        matrix.ncols(),
        [noncommuting_ring(entry.subs({x: 2})) for entry in matrix.list()],
    )
    assert transition_one * transition_two != transition_two * transition_one
    legacy = remainder_forest_pn(
        matrix, moduli, endpoints, 3, kbase=1, indices=indices,
        V=initial, kappa=1,
    )
    def direct_result(modulus, endpoint, nP=3):
        output_ring = PolynomialRing(Integers(modulus), "P")
        current = Matrix(
            output_ring,
            initial.nrows(),
            initial.ncols(),
            [output_ring([
                output_ring(entry.subs({x: 0}))[degree]
                for degree in range(nP)
            ]) for entry in initial.list()],
        )
        for transition in range(1, endpoint):
            step = Matrix(
                output_ring,
                matrix.nrows(),
                matrix.ncols(),
                [output_ring(entry.subs({x: transition})) for entry in matrix.list()],
            )
            product = current * step
            current = Matrix(
                output_ring,
                product.nrows(),
                product.ncols(),
                [output_ring([entry[degree] for degree in range(nP)])
                 for entry in product.list()],
            )
        return current

    direct = {
        index: direct_result(modulus, endpoint)
        for index, modulus, endpoint in zip(indices, moduli, endpoints)
    }
    assert legacy == direct
    native = remainder_forest_pn_native(
        matrix, moduli, endpoints, 3, kbase=1, indices=indices,
        V=initial, kappa=1,
    )
    observed = consumer.read_result(native)
    assert observed["abi_version"] == 1
    assert observed["shape"] == native.shape
    assert observed["endpoints"] == list(endpoints)
    assert [int(value) for value in observed["moduli"]] == moduli
    expected = []
    for index in indices:
        for degree in range(3):
            for row in range(initial.nrows()):
                for column in range(matrix.ncols()):
                    expected.append(int(direct[index][row, column][degree]))
    assert [int(value) for value in observed["coefficients"]] == expected
    assert max(abs(int(value)) for value in observed["coefficients"]) > 2**250

    # Every call creates a separate capsule. A capsule alone keeps its owner
    # and GMP arrays alive after the caller drops the result reference.
    capsule = native.as_capsule()
    second_capsule = native.as_capsule()
    assert second_capsule is not capsule
    later = remainder_forest_pn_native(
        matrix, [97], [5], 3, kbase=1, V=initial, kappa=1
    )
    assert consumer.read_result(later)["shape"] == later.shape
    assert consumer.read_capsule(capsule)["coefficients"] == observed["coefficients"]
    low_precision = remainder_forest_pn_native(
        matrix, [moduli[0]], [endpoints[0]], 1, kbase=1,
        indices=("low",), V=initial, kappa=1,
    )
    low_observed = consumer.read_result(low_precision)
    low_direct = direct_result(moduli[0], endpoints[0], nP=1)
    low_expected = [
        int(low_direct[0, column][0]) for column in range(matrix.ncols())
    ]
    assert low_observed["shape"] == (1, 1, initial.nrows(), matrix.ncols())
    assert [int(value) for value in low_observed["coefficients"]] == low_expected
    del native
    import gc

    gc.collect()
    held = consumer.read_capsule(capsule)
    assert held["coefficients"] == observed["coefficients"]
    assert consumer.read_capsule(second_capsule)["coefficients"] == observed["coefficients"]

    empty = remainder_forest_pn_native(
        matrix, [], [], 2, kbase=1, V=initial
    )
    empty_observed = consumer.read_result(empty)
    assert empty_observed["shape"] == (0, 2, initial.nrows(), matrix.ncols())
    assert empty_observed["endpoints"] == []
    assert empty_observed["moduli"] == []
    assert empty_observed["coefficients"] == []

    modulus_one = remainder_forest_pn_native(
        matrix, [1], [4], 3, kbase=1, V=initial, kappa=1
    )
    modulus_one_observed = consumer.read_result(modulus_one)
    assert set(modulus_one_observed["coefficients"]) == {b"0"}

    import ctypes

    capsule_new = ctypes.pythonapi.PyCapsule_New
    capsule_new.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
    capsule_new.restype = ctypes.py_object
    wrong_name_capsule = capsule_new(
        ctypes.c_void_p(1), b"wrong.name", ctypes.c_void_p()
    )
    try:
        consumer.read_capsule(wrong_name_capsule)
    except ValueError:
        pass
    else:
        raise AssertionError("consumer accepted a capsule with the wrong name")
    print("GMP ABI build tuple: %r" % (observed["gmp"],))
    print("PASS packaged PN GMP ABI compiled consumer; exact large-modulus values and capsule lifetime")


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--isolated-consumer":
        _run_isolated_consumer(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 2:
        build_consumer(sys.argv[1])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
