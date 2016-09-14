all: install

deps: FORCE
		pip install --requirement requirements.txt
	
install: deps
		python setup.py install

test: FORCE
		python setup.py nosetests --attr !long,!huge

long-test: FORCE
		python setup.py nosetests --attr !huge -x

acceptance-huge: FORCE
		python setup.py nosetests

publish: FORCE
		python setup.py sdist upload

gh-pages:
		cd doc; make clean html
		touch doc/_build/html/.nojekyll
		git add doc/
		git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`"
		git subtree split --prefix doc/_build/html -b gh-pages
		git push -f origin gh-pages:gh-pages
		git branch -D gh-pages

clean: FORCE
	rm -rf build/ *.pyc dammit/*.pyc dammit/*.egg-info dammit/*.so dammit/*.c *.egg-info

FORCE:
