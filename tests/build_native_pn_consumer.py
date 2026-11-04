"""Build and exercise a separate consumer using a built pyrforest wheel.

Usage: sage -python tests/build_native_pn_consumer.py path/to/pyrforest.whl
"""

import importlib
import sys
import tempfile
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


def build_consumer(wheel_path):
    wheel_path = Path(wheel_path).resolve()
    source = Path(__file__).with_name("native_pn_consumer.pyx").resolve()
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
        sys.path.insert(0, command.build_lib)
        try:
            return importlib.import_module("native_pn_consumer")
        finally:
            sys.path.remove(command.build_lib)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    consumer = build_consumer(sys.argv[1])

    from sage.all import Integers, Matrix, PolynomialRing, ZZ
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


if __name__ == "__main__":
    main()
