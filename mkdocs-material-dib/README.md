# Material for MkDocs: DIB Lab Theme

A Material Design theme for [MkDocs][1], modified for use by the DIB Lab at UC Davis.

[![Material for MkDocs](docs/assets/images/material.png)][2]

  [1]: http://www.mkdocs.org
  [2]: https://squidfunk.github.io/mkdocs-material/

See [Original Readme](/ORIGINAL_README.md)


## What is MkDocs? What is Material?

[mkdocs](http://www.mkdocs.org/) is a very straightforward 
documentation generator; it is like Sphinx, but for Markdown.

To use mkdocs, you need an `mkdocs.yml` config file
and a `docs/` directory full of Markdown files.
This will generate the site HTML in the `site/` 
directory. (This is all configurable.)

[Material Design](https://material.io/guidelines/material-design/) 
is a set of design principles originating from Google. 
This material theme will give you documentation that 
has a similar look and feel to Google documentation, 
e.g., [Google Cloud ML Engine Documentation](https://cloud.google.com/ml-engine/docs/tensorflow/getting-started-training-prediction).

[Material for Mkdocs](https://squidfunk.github.io/mkdocs-material/)
is an MkDocs theme that applies the principles of Material Design
to generate documentation that looks awesome with minimal effort.


## Quick Start

To get started, you'll follow these steps:

* Create an `mkdocs.yml` config file
* Check out the mkdocs-material-dib theme (optionally: as a submodule)
* (Optional) Set up push to deploy on Github Pages

### MkDocs Config File

Create an `mkdocs.yml` config file:

```
site_name: your-cool-site
site_url: <url-for-your-docs>
repo_name: <repo-owner>/<repo-name>
repo_url: https://github.com/<repo-owner>/<repo-name>
edit_uri: ""

copyright: 'Copyright &copy; 2018 <a href="http://ivory.idyll.org/lab/">Lab for Data Intensive Biology</a> at UC Davis'

# change directory names here
docs_dir: docs
site_dir: site

theme:

  # the next 2 lines are required to use the DIB version
  name: null
  custom_dir: 'mkdocs-material-dib/material'

  # pretty colors! see https://squidfunk.github.io/mkdocs-material/getting-started/#primary-colors
  palette:
    primary: 'blue'
    accent: 'blue'
  
  # fun logos! see https://material.io/icons/
  logo:
    icon: 'dns'

  font:
    text: 'Roboto'
    code: 'Roboto Mono'

# this will add docs/css/custom.css to all your docs
extra_css:
  - css/custom.css


# give a title for each page
nav:
  - Index : 'index.md'
  - About : 'about.md'
  - Special Docs: 'special/docs.md'
  # etc...

# Include the following to enable disqus and
# hypothesis on your mkdocs site:
extra:
  # annotation
  hypothesis: true
  # disqus
  disqus: <insert-name-of-disqus-forum-here>
```

Note: see [mkdocs.yml in the mkdocs-material repo](https://github.com/squidfunk/mkdocs-material/blob/master/mkdocs.yml)
for a more extensive example.

### Install mkdocs

Tags v1.0 and v1.1 of mkdocs-material-dib work with mkdocs <= 0.17.

Tag v2.0 of mkdocs-material-dib works with mkdocs >= 1.0.

We recommend using mkdocs >= 1.0, which can be installed using pip:

```
pip install mkdocs>=1.0
```

### Check Out `mkdocs-material-dib` Theme

Run these commands from wherever you added `mkdocs.yml`
(probably your main repo directory).

Option 1 (no frills): clone a copy of `mkdocs-material-dib`

```
git clone https://github.com/dib-lab/mkdocs-material-dib.git
echo "mkdocs-material-dib" >> .gitignore
echo "site/" >> .gitignore
```

and commit:

```
git add .gitignore
git commit .gitignore -m 'Ignore mkdocs-material-dib'
```

Option 2 (distributable): add `mkdocs-material-dib` as a submodule

```
git submodule add https://github.com/dib-lab/mkdocs-material-dib.git
echo "site/" >> .gitignore
```

and commit:

```
git add .gitmodules .gitignore mkdocs-material-dib
git commit .gitmodules .gitignore mkdocs-material-dib -m 'Add mkdocs-material-dib submodule'
```

**NOTE:** if you have a submodule in a repo, you need to clone it
using the `--recursive` flag to also clone the contents of the 
submodule:

```
git clone --recursive <git-repo-url>
```

If you already have a copy of the repo checked out, but the submodule is empty,
initialize the submodule contents:

```
git submodule update --init
```

### Set Up Push-To-Deploy on Github Pages

The following process will make your final documentation 
available on Github Pages via the following URL:

```
https://<repo-owner>.github.io/<repo-name>
```

When we generate documentation with mkdocs, it puts the 
final static content into a directory called `site/`.

We want to make `site/` into a copy of 
our repo - specifically, the `gh-pages` branch,
which will be a completely separate branch from
the master branch. 

This is a push-to-deploy model: when we push content 
to the `gh-pages` branch, it is deployed to Github Pages.

From the root of this repository, start by removing the site 
directory if it already exists. Now clone a second copy of 
the repo to the `site/` directory:

```
git clone https://github.com/<repo-owner>/<repo-name> site/
cd site/
```

From the site directory, create an orphan branch 
called `gh-pages`. An orphan branch shares no history
with any other branches in the repo. This keeps our 
web content from mixing with our code.

```
git checkout --orphan gh-pages
```

Now remove all the files in the site directory,
and add a "Hello World" page to test that 
Github Pages is working:

```
rm -rf * .gitmodules .gitignore
echo '<h1>Hello World</h1>' > index.html
git add index.html 
git commit index.html -m 'Initial commit of gh-pages branch'
git push origin gh-pages
```

Now go to your repository settings on Github
and enable Github Pages for your repository.
Select the `gh-pages` branch for the content 
location.

Check your Github Pages URL for that hello world page:

```
https://<repo-owner>.github.io/<repo-name>
```

Once you see Hello World, you're ready to upload your docs.
Back in the main directory of your repository,
or wherever `mkdocs.yml` is located, run the command
to build your documentation with mkdocs:

```
mkdocs build
```

This will build all of your documentation 
in the `site/` directory.

To test your documentation out locally, run

```
mkdocs serve
```

and visit `localhost:8000` in your browser.

When you are happy with your documentation,
make a clean copy in preparation for uploading it:

```
rm -rf site/*   # Clean out any existing files
mkdocs build    # Build the docs in site/
cd site/
git add -A .    # Stage every change in the current directory for commit
git commit -a   # Commit all changes
git push origin gh-pages
```

Alternatively, you can use ghp-import:

1. update the docs and commit your changes to master/branch of choice
2. `mkdocs build` to update the docs
3. if you haven’t already, install ghp-import: `conda install -c conda-forge ghp-import`
4. use ghp-import to push the updated to docs to the gh-pages branch  `ghp-import site -p`

Now your documentation should be live!

To re-make your documentation, just re-run the 
block of commands above. 

## Original MkDocs-Material License

Original Github repo: [https://github.com/squidfunk/mkdocs-material](https://github.com/squidfunk/mkdocs-material)

**MIT License**

Copyright (c) 2016-2018 Martin Donath

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
