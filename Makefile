all: install

deps: FORCE
		pip install --requirement requirements.txt
	
install: deps
		python setup.py install

test: FORCE
		py.test -m "not long and not huge"

ci-test: FORCE
		py.test -m "not long and not huge and not requires_databases"

long-test: FORCE
		py.test -m "not huge"

acceptance-huge: FORCE
		py.test 

publish: FORCE
		python setup.py sdist upload

doc: FORCE
		sphinx-apidoc  -o doc/ dammit/ -f
		cd doc; make clean html

banner:
		toilet --gay -f big dammit! > dammit/banner.txt

gh-pages: doc
		touch doc/_build/html/.nojekyll
		git add doc/
		git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`"
		git subtree split --prefix doc/_build/html -b gh-pages
		git push -f origin gh-pages:gh-pages
		git branch -D gh-pages

clean: FORCE
	rm -rf build/ *.pyc dammit/*.pyc dammit/*.egg-info dammit/*.so dammit/*.c *.egg-info

FORCE:
