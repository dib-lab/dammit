#!/usr/bin/env python

import argparse
from conda_build import api
from conda_build.environ import create_env
import conda.cli


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('recipe_path')
    args = parser.parse_args()

    metadata = api.render(args.recipe_path)[0][0]
    deps = metadata.get_value('requirements/run') \
           + metadata.get_value('test/requires')

    conda.cli.main('conda', 'install',  '-y', ' '.join(deps))

if __name__ == '__main__':
    main()
