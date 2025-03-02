from packaging import version

from quantmsrescore.openms import OpenMSHelper


def test_get_pyopenms_version():
    oms_version = OpenMSHelper.get_pyopenms_version()
    installed_version = version.parse(oms_version)
    assert installed_version >= version.parse("2.8.0")
