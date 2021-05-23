#!/usr/bin/env python


def ductile_test(ratio):
    """Test for ductility."""
    if ratio > 1.75:
        return "ductile"
    else:
        return "brittle"
