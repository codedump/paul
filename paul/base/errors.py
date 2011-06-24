
"""Provides exceptions for Paul Loader."""

#
# Wrong file version (which means that file format was ok)
#
class VersionError (Exception):
    def __init__ (self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

#
# Wrong file format
#
class FormatError (Exception):
    def __init__ (self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

#
# Feature not implemented
#
class NotImplementedError(Exception):
    def __init__ (self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
