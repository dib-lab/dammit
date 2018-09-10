#!/usr/bin/env python

import argparse
import os
import sys

from conda_build import api
from conda_build.environ import get_install_actions
from conda_build.index import get_build_index
from conda_build.conda_interface import display_actions, execute_actions
import conda.cli


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('recipe_path')
    parser.add_argument('--prefix', default=os.environ['CONDA_PREFIX'])
    parser.add_argument('--env', default=os.environ['CONDA_DEFAULT_ENV'])
    args = parser.parse_args()

    metadata = api.render(args.recipe_path)[0][0]
    config = metadata.config
    deps = metadata.get_value('requirements/build') \
           + metadata.get_value('requirements/run') \
           + metadata.get_value('test/requires')
    deps = list(set(deps))
    
    actions = get_install_actions(args.prefix, tuple(deps), args.env,
                              subdir=metadata.config.subdir,
                              verbose=config.verbose,
                              debug=config.debug,
                              locking=config.locking,
                              bldpkgs_dirs=tuple(config.bldpkgs_dirs),
                              timeout=config.timeout,
                              disable_pip=config.disable_pip,
                              max_env_retry=config.max_env_retry,
                              output_folder=config.output_folder,
                              channel_urls=tuple(config.channel_urls))

    index, index_ts = get_build_index(subdir=metadata.config.subdir,
                                    bldpkgs_dir=config.bldpkgs_dir,
                                    output_folder=config.output_folder,
                                    channel_urls=config.channel_urls,
                                    debug=config.debug,
                                    verbose=config.verbose,
                                    locking=config.locking,
                                    timeout=config.timeout)

    display_actions(actions, index)
    execute_actions(actions, index, verbose=config.debug)


if __name__ == '__main__':
    main()
