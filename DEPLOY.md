# How to deplot to PyPI

1. Update the version number in `pyproject.toml` to match the `__version__` number in the `src/` directory.
2. Clean up the previous builds, make the new build and upload. You may need to first run `pip install --upgrade build twine`.

```python
rm -rf dist/ build/ *.egg-info src/*.egg-info*    # Clean up old builds first
python -m build                   # Build the new version
python -m twine upload dist/*     # Upload to PyPI
```

3. Double check on PyPI that the new version uploaded correctly.
