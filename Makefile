all:
		install

deps: FORCE
		pip install --requirement requirements.txt
	
install: deps
		python setup.py install

test: FORCE
		python setup.py test

publish: FORCE
		python setup.py sdist upload

gh-pages:
		cd doc; make clean html
		touch doc/_build/html/.nojekyll
		git subtree push --prefix doc/_build/html origin gh-pages
		git add doc/
		git ci -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout master

clean: FORCE
	rm -rf build/ *.pyc dammit/*.pyc dammit/*.egg-info dammit/*.so dammit/*.c *.egg-info

FORCE:
