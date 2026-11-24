def pytest_configure(config):
    config.addinivalue_line("markers", "slow: re-solves models (a few minutes); deselect with -m 'not slow'")
