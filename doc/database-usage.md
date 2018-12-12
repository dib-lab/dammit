---
title: 'Database Usage'
---

dammit handles databases under the `dammit databases` subcommand. By
default, dammit looks for databases in
`$HOME/.dammit/databases` and will install them there if
missing. If you have some of the databases already, you can inform
dammit with the `--database-dir` flag.

To check for databases in the default location:

```
dammit databases
```

To check for them in a custom location, you can either use the
`--database-dir` flag:

```
dammit databases --database-dir /path/to/databases
```

or, you can set the `DAMMIT_DB_DIR` environment variable.
The flag will supersede this variable, falling back to the default if
neither is set. For example:

```
export DAMMIT_DB_DIR=/path/to/databases
```

This can also be added to your `$HOME/.bashrc` file to make
it persistent.

To download and install them into the default directory:

```
dammit databases --install
```

For info on the specific databases used in dammit, see [About Databases](database-about.md).

For advanced installation and usage instructions, check out the
[Advanced Database Handling](database-advanced.md) section.
