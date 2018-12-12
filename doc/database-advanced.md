---
title: 'Advanced Database Handling'
---

Several of these databases are quite large. Understandably, you probably
don't want to download or prepare them again if you already have. There
are a few scenarios you might run in to.

1.  You already have the databases, and they're all in one place and
    properly named.

    > Excellent! This is the easiest. You can make use of dammit's
    > `--database-dir` flag to tell it where to look. When running with
    > `--install`, it will find the existing files and prep them if
    > necessary.:
    >
    >     dammit databases --database-dir <my_database_dir> --install

2.  Same as above, but they have different names.

    > dammit expects the databases to be "properly" named -- that is,
    > named the same as their original forms. If your databases aren\'t
    > named the same, you'll need to fix them. But that's okay! We can
    > just soft link them. Let's say you have Pfam-A already, but for
    > some reason its named `all-the-models.hmm`. You can
    > link them to the proper name like so:
    >
    >     cd <my_database_dir>
    >     ln -s all-the-models.hmm Pfam-A.hmm
    >
    > If you already formatted it with `hmmpress`, you can
    > avoid repeating that step as well:
    >
    >     ln -s all-the-models.hmm.h3f Pfam-A.hmm.h3f
    >     ln -s all-the-models.hmm.h3i Pfam-A.hmm.h3i
    >     ln -s all-the-models.hmm.h3m Pfam-A.hmm.h3m
    >     ln -s all-the-models.hmm.h3p Pfam-A.hmm.h3p
    >
    > For a complete listing of the expected names, just run the
    > `databases` command:
    >
    >     dammit databases

3.  You have the databases, but they're scattered to the virtual winds.

    > The fix here is similar to the above. This time, however, we'll
    > soft link all the databases to one location. If you've run
    > `dammit databases`, a new directory will have been created at
    > `\$HOME/.dammit/databases`. This is where they are
    > stored by default, so we might as well use it! For example:
    >
    >     cd $HOME/.dammit/databases
    >     ln -s /path/to/all-the-models.hmm Pfam-A.hmm
    >
    > And repeat for all the databases. Now, in the future, you will be
    > able to run dammit without the `--database-dir` flag.

Alternatively, if this all seems like too much of a hassle and you have
lots of hard drive space, you can just say "to hell with it!" and
reinstall everything with:

```
dammit databases --install
```

