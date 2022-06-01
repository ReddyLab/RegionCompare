.PHONY: clean install

PYTHON_FILES=$(wildcard RegionCompare/*/*.py)

clean:
	-rm -rf RegionCompare.egg-info build dist

dist: $(PYTHON_FILES)
	if [ -z $$(pip list | grep -e "^build\s") ]; then pip install build; fi
	python -m build

install: dist
	pip install --force-reinstall dist/*.whl
