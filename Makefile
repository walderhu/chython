test:
	python3 testing.py

gcov:
	python3 -m coverage run -m pytest unit.py
	python3 -m coverage report
	python3 -m coverage html

clean: 
	rm -rf __pycache__ .pytest_cache htmlcov .coverage output_files
	rm -rf gcov/.pytest_cache  gcov/htmlcov gcov/.coverage gcov/output_files 

project:
	pip install CachedMethods numpy lazy_object_proxy importlib_resources lxml

install_sympy:
	sudo apt install python3-sympy

install:
	pip install coverage
	pip install pytest

push:
	git add * && git commit && git push origin $(shell git rev-parse --abbrev-ref HEAD)
	pip install pyppeteer