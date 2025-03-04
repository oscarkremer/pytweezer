.DEFAULT_GOAL = check
PYTHON_INTERPRETER = python3
PROJECT_NAME = pytweezer
################################################################################
# COMMANDS                                                                     #
################################################################################


install: ##Installation method. Creates folder named as 'dirs' and conda environment with dependences. \n\
		Run verification and configuration steps for requisites.
install: test_environment dirs
	@echo "---> Running setup ..."
	@conda env create -q -f environment.yml --name $(PROJECT_NAME) > /dev/null
	@echo "---> To complete setup please run: conda activate $(PROJECT_NAME)"


test_environment: ##Run verification step for requisites
	@echo "---> Verifying if requisites are already installed"
	@$(PYTHON_INTERPRETER) test_environment.py


update: ##Update method to change dependences.
	@echo "---> Updating dependencies"
	@conda env update -q -f environment.yml

set:
	@python setup.py install


dirs:	## Make command to create folder for results dataframes.
	@echo "---> Creating data folder for results"
	@mkdir -p data/results
	@mkdir -p data/logs
	@echo "----> Done"

test:
	@coverage run -m pytest tests/ -v

help:	## Help method to list all available commands.
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'
 

clean:	## Method for removing cached and .pyc or .pyo files.
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

upload:
	python -m twine upload dist/*
