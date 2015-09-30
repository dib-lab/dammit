all:
		python setup.py install

test:
		python setup.py test

publish:
		python setup.py sdist upload

clean:
	rm -rf build/ *.pyc dammit/*.pyc dammit/*.egg-info dammit/*.so dammit/*.c *.egg-info
