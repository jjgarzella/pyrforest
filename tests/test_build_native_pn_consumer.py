"""Negative controls for the packaged PN consumer's producer provenance."""

import importlib.util
import tempfile
from pathlib import Path


HELPER_PATH = Path(__file__).with_name("build_native_pn_consumer.py")
SPEC = importlib.util.spec_from_file_location("build_native_pn_consumer", HELPER_PATH)
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)


def test_wheel_package_and_compiled_extension_paths_are_accepted():
    with tempfile.TemporaryDirectory() as temporary:
        package = Path(temporary) / "wheel" / "pyrforest"
        package.mkdir(parents=True)
        package_file = package / "__init__.py"
        extension_file = package / "rforest.cpython-312-darwin.so"
        package_file.write_text("# wheel package\n")
        extension_file.write_bytes(b"test extension")
        assert helper._validated_producer_paths(
            package_file, extension_file, package
        ) == (package_file.resolve(), extension_file.resolve())


def test_ambient_package_outside_wheel_is_rejected():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        wheel_package = root / "wheel" / "pyrforest"
        ambient_package = root / "ambient" / "pyrforest"
        wheel_package.mkdir(parents=True)
        ambient_package.mkdir(parents=True)
        (wheel_package / "__init__.py").write_text("# wheel package\n")
        ambient_init = ambient_package / "__init__.py"
        ambient_extension = ambient_package / "rforest.cpython-312-darwin.so"
        ambient_init.write_text("# conflicting producer\n")
        ambient_extension.write_bytes(b"conflicting extension")
        try:
            helper._validated_producer_paths(
                ambient_init, ambient_extension, wheel_package
            )
        except AssertionError as error:
            assert "outside the expected wheel directory" in str(error)
        else:
            raise AssertionError("ambient producer paths were accepted")


def test_source_file_cannot_stand_in_for_compiled_extension():
    with tempfile.TemporaryDirectory() as temporary:
        package = Path(temporary) / "wheel" / "pyrforest"
        package.mkdir(parents=True)
        package_file = package / "__init__.py"
        source_extension = package / "rforest.py"
        package_file.write_text("# wheel package\n")
        source_extension.write_text("# wrong producer module\n")
        try:
            helper._validated_producer_paths(
                package_file, source_extension, package
            )
        except AssertionError as error:
            assert "compiled extension" in str(error)
        else:
            raise AssertionError(
                "a Python source file was accepted as the producer extension"
            )


def test_missing_producer_extension_is_rejected():
    with tempfile.TemporaryDirectory() as temporary:
        package = Path(temporary) / "wheel" / "pyrforest"
        package.mkdir(parents=True)
        package_file = package / "__init__.py"
        package_file.write_text("# wheel package\n")
        try:
            helper._validated_producer_paths(package_file, None, package)
        except AssertionError as error:
            assert "was not loaded" in str(error)
        else:
            raise AssertionError("a missing producer extension was accepted")


def main():
    test_wheel_package_and_compiled_extension_paths_are_accepted()
    test_ambient_package_outside_wheel_is_rejected()
    test_source_file_cannot_stand_in_for_compiled_extension()
    test_missing_producer_extension_is_rejected()
    print("PASS PN consumer producer-provenance positive and negative controls")


if __name__ == "__main__":
    main()
