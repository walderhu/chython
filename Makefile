FILE = main.py
TEST_FILE = unit.py
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)

push:
	git add * && git commit && git push origin $(BRANCH)

install_sympy:
	sudo apt install python3-sympy

install:
	pip install coverage
	pip install pytest

test:
	python3 $(TEST_FILE)

gcov:
	python3 -m coverage run -m pytest $(TEST_FILE)
	python3 -m coverage report
	python3 -m coverage html
#	open htmlcov/index.html

clean: 
	rm -rf __pycache__ .pytest_cache htmlcov .coverage

clang:
	@autopep8 --in-place --aggressive --aggressive $(FILE)
	@autopep8 --in-place --aggressive --aggressive $(TEST_FILE)
