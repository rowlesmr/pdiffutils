import pdiffutils as pd


def test_val_err_str():
    assert pd.val_err_str(12.345, 1.23) == "12.3(12)"
    assert pd.val_err_str(12.345, 0.0123) == "12.345(12)"
    assert pd.val_err_str(12345, 123) == "12340(120)"








